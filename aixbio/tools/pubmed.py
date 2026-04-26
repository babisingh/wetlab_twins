from __future__ import annotations

import logging
import time
import xml.etree.ElementTree as ET

import httpx

logger = logging.getLogger(__name__)

_ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_TIMEOUT = 12.0


def search_expression_literature(
    protein_name: str,
    uniprot_id: str,
    max_results: int = 5,
) -> list[dict]:
    """Query PubMed for E. coli expression papers about this protein.

    Returns a list of dicts with keys: pmid, title, abstract, year.
    Returns [] on any network or parse error (non-fatal).
    """
    query = (
        f'"{protein_name}"[Title/Abstract] AND "expression"[Title/Abstract] '
        f'AND ("E. coli"[Title/Abstract] OR "Escherichia coli"[Title/Abstract])'
    )

    try:
        pmids = _esearch(query, max_results)
        if not pmids:
            logger.info("PubMed: no results for %s (%s)", protein_name, uniprot_id)
            return []

        time.sleep(0.4)  # polite delay between Entrez requests
        articles = _efetch(pmids)
        logger.info("PubMed: retrieved %d articles for %s", len(articles), protein_name)
        return articles

    except Exception as exc:
        logger.warning("PubMed fetch failed for %s: %s", protein_name, exc)
        return []


def _esearch(query: str, max_results: int) -> list[str]:
    resp = httpx.get(
        f"{_ENTREZ_BASE}/esearch.fcgi",
        params={
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "retmode": "json",
            "sort": "relevance",
        },
        timeout=_TIMEOUT,
    )
    resp.raise_for_status()
    return resp.json().get("esearchresult", {}).get("idlist", [])


def _efetch(pmids: list[str]) -> list[dict]:
    resp = httpx.get(
        f"{_ENTREZ_BASE}/efetch.fcgi",
        params={
            "db": "pubmed",
            "id": ",".join(pmids),
            "rettype": "abstract",
            "retmode": "xml",
        },
        timeout=_TIMEOUT,
    )
    resp.raise_for_status()
    return _parse_xml(resp.text)


def _parse_xml(xml_text: str) -> list[dict]:
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []

    articles = []
    for article in root.findall(".//PubmedArticle"):
        pmid = _text(article, ".//PMID")
        title = _text(article, ".//ArticleTitle")
        abstract = _text(article, ".//AbstractText")
        year = _text(article, ".//PubDate/Year") or _text(article, ".//PubDate/MedlineDate")[:4]

        if pmid and (title or abstract):
            articles.append({
                "pmid": pmid,
                "title": title,
                "abstract": abstract[:700],  # keep tokens manageable
                "year": year,
            })
    return articles


def _text(element, xpath: str) -> str:
    el = element.find(xpath)
    if el is None:
        return ""
    return (el.text or "").strip()

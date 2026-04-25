from __future__ import annotations

import os

DEFAULT_HOST = "E. coli K12"
DEFAULT_AVOID_SITES = ("BamHI", "XhoI", "EcoRI")
DEFAULT_TAG_TYPE = "6xHis"
DEFAULT_PROTEASE_SITE = "Enterokinase"
DEFAULT_VECTOR = "pET-28a(+)"
DEFAULT_CLONING_SITES = ("BamHI", "XhoI")
DEFAULT_MAX_REMEDIATION_ATTEMPTS = 3
LLM_MODEL = os.getenv("LLM_MODEL", "anthropic/claude-sonnet-4-6")
OPENROUTER_BASE_URL = "https://openrouter.ai/api/v1"
OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY", "")

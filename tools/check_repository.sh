#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

status=0
matlab_count="$(find code -type f -name '*.m' | wc -l)"
echo "MATLAB files: $matlab_count"

if [[ "$matlab_count" -ne 194 ]]; then
    echo "ERROR: expected 194 MATLAB files."
    status=1
fi

if rg -n -i --glob '*.m' --glob '*.md' \
    '(BEGIN [A-Z ]*PRIVATE KEY|api[_-]?key[[:space:]]*=|password[[:space:]]*=|github_pat_[A-Za-z0-9_]+|ghp_[A-Za-z0-9]+)' .; then
    echo "ERROR: possible credential material detected."
    status=1
else
    echo "Credential scan: clear"
fi

path_count="$(rg -l --glob '*.m' '(/home/|/Users/|[A-Za-z]:\\)' code | wc -l || true)"
echo "Historical files with machine-specific paths: $path_count"

bad_names="$(find code -type f -name '*.m' ! -regex '.*/\(Photon_[0-9][0-9]\([_][A-Za-z0-9]*\)?\|Tartrazine_[0-9][0-9]\)\.m' -print)"
if [[ -n "$bad_names" ]]; then
    echo "ERROR: unexpected MATLAB filenames:"
    echo "$bad_names"
    status=1
fi

if git check-ignore -q config/project_config.example.m; then
    echo "ERROR: example configuration is unexpectedly ignored."
    status=1
fi

if [[ -e config/project_config.m ]] && ! git check-ignore -q config/project_config.m; then
    echo "ERROR: local configuration is not ignored."
    status=1
fi

exit "$status"

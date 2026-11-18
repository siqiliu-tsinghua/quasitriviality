# Security policy

This is a local research library, not a sandbox or a network service.
The public calculation API takes constructed SymPy expressions; it does not
parse user-supplied Python or SymPy expression strings. Do not pass untrusted
Python objects or callbacks to the API, or parse untrusted text using an
evaluator upstream. The example tools accept local JSON coefficient archives;
use trusted archives. Symbolic inputs or archives can consume substantial
CPU and memory. If you build a service around the library, validate inputs
and impose independent size, time, and memory limits.

Report security findings privately via **Security → Advisories → Report a
vulnerability** on the GitHub repository, if available. Otherwise ask a
maintainer for a private reporting channel without posting the details in a
public issue.

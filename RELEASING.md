# Release checklist

This project is a finite-order research prototype. Adding research examples
is not a prerequisite to publishing the source.

Before making a GitHub repository public:

1. Review the files and all branches or tags to be pushed, including commit
   author and committer information. A clean working tree does not remove
   information from existing Git history. Check the intended public
   attribution; GitHub has a
   [commit email privacy guide](https://docs.github.com/en/account-and-profile/how-tos/email-preferences/setting-your-commit-email-address).
2. Keep the MIT license and reference attributions. No reference PDFs, Lean
   source, or Mathematica source are bundled. The benchmark archives contain
   computed results; literature formulas are cited where compared.
3. Run `python -m unittest discover -v` from the source checkout, and run
   the README quick start. For a fresh high-order verification use
   [Toda's independent reproduction commands](examples/toda/HIGH_ORDER.md)
   and the [loop-equation check](examples/toda/LOOP_EQUATION.md).
4. If distributing a Python package, build the wheel and source distribution
   from a clean checkout, inspect their contents, and verify that they
   include `LICENSE`. The wheel installs the library; use the checkout or
   source distribution for the benchmark runners and data. Check the
   advisories for the build and runtime dependencies selected for the release.
5. Preview the README using [GitHub's Markdown math syntax](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions).
   If using GitHub for issue reporting, enable private vulnerability reporting
   so the route in [SECURITY.md](SECURITY.md) is available.

Publishing the source on GitHub is separate from uploading packages to PyPI.

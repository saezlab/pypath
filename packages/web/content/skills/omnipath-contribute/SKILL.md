---
name: omnipath-contribute
description: Contribute a new biological data resource to OmniPath by adding an inputs_v2 module to saezlab/pypath and opening a pull request.
---

# Contribute to OmniPath

Add the user's chosen data resource to [saezlab/pypath](https://github.com/saezlab/pypath)
as a new `pypath/inputs_v2` module, validate it, and open an upstream pull request.
If the resource is unspecified, ask for its name or source URL.

Use an existing pypath checkout or clone the repository. In the OmniPath development
workspace, the checkout is commonly `../omnipath-build/pypath`. Follow its `AGENTS.md`
and verify the appropriate upstream base branch for `inputs_v2` contributions.

Before implementing, read:

- `pypath/inputs_v2/README.md` for the declarative input framework and contribution workflow.
- `pypath/inputs_v2/BIOLINK_MODELING.md` for biological modeling conventions.
- Existing `inputs_v2` modules and their tests that resemble the new source's format
  and biological content. Use these as implementation examples.

Treat those documents and the current code as the source of truth. Inspect the
resource's official documentation, data format and data license, then implement the
module and any necessary parser or resource metadata. Follow the documented validation
steps and add representative tests using the repository's current test setup.

Work on a dedicated branch, preserving unrelated changes. Commit the contribution
and open a pull request against `saezlab/pypath`, using a fork if needed. Describe the
imported datasets, source and license, checks actually run, and any known limitations.
Use a draft PR if validation remains incomplete. Return the PR URL; if access prevents
submission, leave the commit and PR description ready and explain what remains.

Respect a request to stop at a local patch. Do not merge the PR or deploy OmniPath
as part of adding an input module.

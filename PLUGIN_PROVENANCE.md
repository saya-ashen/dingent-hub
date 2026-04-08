# Plugin Provenance And Risk Notes

This document records known upstream/source information for plugins imported from `ENA/` and `Public_plugins/`, and highlights operational/security risks that should be reviewed before production use.

## Newly Added Plugins

### Remote HTTP MCP Endpoints

- `plugins/ENA/plugin.toml`: remote endpoint `https://nucleotide-archive-mcp.fastmcp.app/mcp`
- `plugins/Omnipath/plugin.toml`: remote endpoint `https://explore.omnipathdb.org/api/mcp`
- `plugins/ClinicalTrials/plugin.toml`: remote endpoint `https://clinicaltrials.caseyjhand.com/mcp`
- `plugins/STRING/plugin.toml`: remote endpoint `https://mcp.string-db.org/`

### Local/Stdio Plugins With Upstream Clues

- `plugins/pubchem/README.md`: references `https://github.com/JackKuo666/PubChem-MCP-Server.git`
- `plugins/uniprot/README.md`: references `https://github.com/TakumiY235/uniprot-mcp-server.git`
- `plugins/cellmarker/README.md`: references docs `https://docs.scmcphub.org` and CellMarker citation DOI `https://doi.org/10.1093/nar/gkac947`
- `plugins/chEMBL/README.md`: contains placeholder clone URL (`https://github.com/yourusername/ChEMBL-MCP-Server.git`), so exact upstream repo is not confirmed

## Risk Checklist

- Remote MCP trust boundary: HTTP plugins send query content to third-party services. Validate data policy before enabling in sensitive environments.
- Availability and change risk: remote services may change API behavior, throttle requests, or go offline.
- Upstream integrity risk: for plugins without explicit, verifiable origin metadata, treat code as unverified until repository and maintainer identity are confirmed.
- License compliance: not all imported plugins clearly state license in plugin metadata. Verify license terms before redistribution/commercial use.
- Dependency risk: stdio plugins rely on third-party dependencies (for example, external API clients). Run dependency and vulnerability scans before release.

## Recommended Controls

- Pin plugin versions and dependency lockfiles where possible.
- Restrict network egress for plugin runtime to required domains only.
- Add periodic health checks for remote MCP endpoints.
- Record owner, upstream repository URL, and license for each plugin in a managed inventory.

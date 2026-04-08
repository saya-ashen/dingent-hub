# dingent-hub

Plugin repository for dingent-compatible MCP plugins.

## Plugin Source And Risk Notes

See `PLUGIN_PROVENANCE.md` for upstream information and risk notes for third-party plugins (especially remote HTTP MCP endpoints and externally sourced plugin implementations).

## Newly added plugins

The following plugins were added from `ENA/` and `Public_plugins/` into `plugins/`:

- `plugins/ENA`
- `plugins/chEMBL`
- `plugins/Omnipath`
- `plugins/ClinicalTrials`
- `plugins/pubchem`
- `plugins/uniprot`
- `plugins/cellmarker`
- `plugins/STRING`

## Source and risk notes

To reduce supply-chain and runtime risks, keep the following source/provenance notes in mind.

### Remote MCP endpoints

These plugins call third-party remote MCP services directly:

- `plugins/ENA/plugin.toml` -> `https://nucleotide-archive-mcp.fastmcp.app/mcp`
- `plugins/Omnipath/plugin.toml` -> `https://explore.omnipathdb.org/api/mcp`
- `plugins/ClinicalTrials/plugin.toml` -> `https://clinicaltrials.caseyjhand.com/mcp`
- `plugins/STRING/plugin.toml` -> `https://mcp.string-db.org/`

Risks and controls:

- Availability depends on external service uptime.
- Response content may change without notice.
- Requests may include user query data; evaluate privacy/compliance before production use.
- Prefer allow-listing trusted domains and adding timeout/retry limits in clients.

### Public code provenance (from plugin READMEs)

- `plugins/pubchem/README.md` references: `https://github.com/JackKuo666/PubChem-MCP-Server.git`
- `plugins/uniprot/README.md` references: `https://github.com/TakumiY235/uniprot-mcp-server.git`
- `plugins/cellmarker/README.md` references docs: `https://docs.scmcphub.org`
- `plugins/chEMBL/README.md` includes a placeholder clone URL (`yourusername`), so upstream origin is not uniquely identified.

Recommended follow-up:

- Pin upstream commit/tag when origin is known.
- Verify license and maintenance status for each imported plugin.
- Track local modifications separately from upstream mirrors.

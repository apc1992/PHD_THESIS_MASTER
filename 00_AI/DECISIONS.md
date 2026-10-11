# ARCHITECTURAL & SCIENTIFIC DECISIONS LOG

| Date | Topic | Decision Made | Rationale | Impacted Files |
|------|-------|---------------|-----------|----------------|
| 2026-09-18 | Initial AI workspace | Created central `00_AI` control layer and project-specific `AI_WORK` areas. | Keep project rules, decisions, evidence tracking, and AI outputs separate from source research material. | All projects |
| 2026-09-21 | Thesis scope | Thesis has 3 scientific parts: current antigen-discovery approaches; pan-cancer platform application; mimotope optimization. | Thesis narrative is broader than either publication and must connect discovery to optimization. | `PROJECT_MAP.md`, `THESIS_CHAPTER_MAP.md` |
| 2026-09-21 | Priority update | Pan-cancer work is active only to the minimum level needed for the thesis; Publication 2 manuscript is deferred. | Parts I–II depend on preliminary pan-cancer analysis, while manuscript completion is lower priority. | `PRIORITIES.md`, `03_PANCAN` |
| 2026-09-21 | Mimotope status | Computational review/remake is active while final immunogenicity experiments are pending. | The thesis Part III can be developed before all experiments are complete, provided pending evidence is clearly labeled. | `01_MIMOTOPES`, thesis Part III |
| 2026-09-21 | Data organization | Use a sample × analysis workbook as the central pan-cancer inventory. | A matrix exposes analysis coverage/gaps before new work is proposed. | `03_PANCAN/cohort` |
| 2026-09-21 | Execution environment | Local Mac is the default for writing, literature, lightweight coding and figures; HPC is reserved for heavy computation as needed. | Reduces infrastructure complexity during the thesis deadline period. | All computational work |

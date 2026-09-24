# Capture templates

- `shot_template.json`: copy to a per-shot file, replace every `REPLACE-...`
  value, and either fill inline `samples` or delete `samples` and supply an
  external CSV via `--series`.
- `series_template.csv`: header
  `time_s,beverage_mass_g,pressure_Pa,pressure_location,temperature_K,temperature_location`.
  Leave unknown optionals empty. Do not invent pressure location; a boiler
  gauge is not puck pressure.

Keep `source` as `measured` for real shots, `synthetic` for tests. Never put
exact machine serials in shared files. Record `time_zero` as
`pump_activation`, `first_liquid`, or `already_wet`.
See `SHOT_DATA.md` for validation rules and the normal brewing workflow.

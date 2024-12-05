---
title: Database Issues
area: deployment
box_type: tip
layout: faq
contributors: [natefoo]
---

For slow queries, start with `EXPLAIN ANALYZE`
  - [Example of the "jobs ready to run" query](https://gist.github.com/natefoo/da46bea8136ba67d65f616c86a27c454)

However it can be useful to dig into the queries with the [Postgres EXPLAIN Visualizer (PEV)](https://tatiyants.com/pev/#/plans) to get a more visual and clear representation. (Try it with this [demo data](https://gist.github.com/hexylena/467c7726b5ab9e18c47080893dbc072e))

You can set some options in the Galaxy configuration or database that will help debugging this:

- `database_engine_option_echo` (but warning, extremely verbose)
- `slow_query_log_threshold` logs to Galaxy log file
- `sentry_sloreq_threshold` if using Sentry

Additionally check that your database is running `VACUUM` regularly enough and look at `VACUUM ANALYZE`

There are some `gxadmin query pg-*` commands which can help you monitor and track this information.

Lastly, check your database settings! It might not have enough resources allocated. Check [PGTune](https://pgtune.leopard.in.ua/) for some suggestions of optimised parameters.

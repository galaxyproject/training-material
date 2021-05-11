---
layout: tutorial_hands_on

title: "Data source integration"
questions:
 - How can I develop extensions to Galaxy data model?
 - How can I implement new API functionality within Galaxy?
 - How can I extend the Galaxy user interface with VueJS components?
objectives:
time_estimation: "180M"
contributors:
 - jmchilton
key_points:
 - Galaxy database interactions are mitigated via SQL Alchemy code in lib/galaxy/model.
 - Galaxy API endpoints are implemented in lib/galaxy/webapps/galaxy, but generally defer to application logic in lib/galaxy/managers. 
 - Galaxy client code should do its best to separate API interaction logic from display components.
---

# Contributed to Galaxy Core 

This tutorial walks you through an extension to Galaxy and how to contribute back to the core project.


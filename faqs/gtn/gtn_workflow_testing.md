---
title: Adding workflow tests with Planemo
area: contributors
layout: faq
box_type: tip
contributors: [hexylena]
---

1. Find a tutorial that you're interested in, that **doesn't currently have tests.**

   This tutorial has a workflow (`.ga`) and a test, notice the `-test.yml` that has the same name as the workflow `.ga` file.

   ```
   machinelearning/workflows/machine_learning.ga
   machinelearning/workflows/machine_learning-test.yml
   ```

   You want to find tutorials without the `-test.yml` file. The workflow file might also be missing.

2. Check if it has a workflow (if it does, skip to step 5.)
3. Follow the tutorial
4. Extract a workflow from the history
5. Run that workflow in a new history to test
6. Obtain the workflow invocation ID, and your API key (User → Preferences → Manage API Key)

   ![screenshot of the workflow invocation page. The user drop down shows where to find this page, and a red box circles a field named "Invocation ID"]({% link faqs/gtn/images/invocation.png %})

7. Install the latest version of `planemo`

   ```
   # In a virtualenv
   pip install planemo
   ```

8. Run the command to initialise a workflow test from the `workflows/` subdirectory - if it doesn't exist, you might need to create it first.

   ```
   planemo workflow_test_init --from_invocation <INVOCATION ID> --galaxy_url <GALAXY SERVER URL> --galaxy_user_key" <GALAXY API KEY>
   ```

   This will produce a folder of files, for example from a testing workflow:

   ```
   $ tree
   .
   ├── test-data
   │   ├── input dataset(s).shapefile.shp
   │   └── shapefile.shp
   ├── testing-openlayer.ga
   └── testing-openlayer-tests.yml
   ```

9. You will need to check the `-tests.yml` file, it has some automatically generated comparisons. Namely it tests that output data matches the test-data exactly, however, you might want to replace that with assertions that check for e.g. correct file size, or specific text content you expect to see.
10. Contribute all of those files to the GTN in a PR.

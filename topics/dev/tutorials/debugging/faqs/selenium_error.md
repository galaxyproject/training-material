---
title: When I try to run a Selenium test, I get an error
box_type: question
layout: faq
contributors: [assuntad23]
---

If you get the following error:

```bash
selenium.common.exceptions.SessionNotCreatedException (...This version of ChromeDriver only supports Chrome version...)
```

Make sure that (a) the version of your ChromeDriver is the same as the version of Chrome:

```bash
$ chromedriver --version
$ chrome --version
```

If they are not the same:
- download the appropriate version of [ChromeDriver](https://chromedriver.chromium.org/downloads).
- unzip the file
- move the chromedriver file into the appropriate location.
  - On Linux, that could be `/usr/bin`, `$HOME/.local/bin`, etc.
  - Use the `which` command to check the location: `$ which chromedriver`
- Make sure the permissions are correct (755).


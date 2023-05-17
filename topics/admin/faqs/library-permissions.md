---
title: Library Permission Issues
area: data-libraries
box_type: tip
layout: faq
contributors: [hexylena]
---

When running `setup-data-libraries` it imports the library with the permissions of the admin user, rather locked down to the account that handled the importing.

Due to how data libraries have been implemented, it isn't sufficient to share the folder with another user, instead you must also share individual items within this folder. This is an unfortunate issue with Galaxy that we hope to fix someday.

Until then, we can recommend you install the latest version of Ephemeris which includes the `set-library-permissions` command which let's you recursively correct the permissions on a data library. Simply run:

```
set-library-permissions -g https://galaxy.example.com -a $API_KEY LIBRARY --roles ROLES role1,role2,role3
```

Where `LIBRARY` is the id of the library you wish to correct.

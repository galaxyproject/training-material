---
title: Upload many files (>10) via FTP
area: data upload
box_type: tip
layout: faq
contributors: [bebatut, shiltemann]
---

Some Galaxies offer FTP upload for very large datasets.

1. Check that your Galaxy supports FTP upload and look up the FTP settings.
   - **Galaxy EU** (usegalaxy.eu): [FTP settings](https://usegalaxy-eu.github.io/ftp/)
   - **Galaxy Australia** (usegalaxy.org.au): [FTP settings](https://usegalaxy-au.github.io/posts/2019/03/18/new-ftp-upload-url/plain.html)
   - **Galaxy Main** (usegalaxy.org): [FTP no longer supported, command-line upload available](https://help.galaxyproject.org/t/the-usegalaxy-org-ftp-service-will-be-decommissioned-on-august-12-2022/8318)
   - **Other Galaxies**: check the homepage of your Galaxy or contact one of the Galaxy administrators.

2. Make sure to have an FTP client installed

    There are many options. We can recommend [FileZilla](https://filezilla-project.org/), a free FTP client that is available on Windows, MacOS, and Linux.

3. Establish FTP connection to the Galaxy server
    1. Provide the Galaxy server's FTP server name (e.g. `ftp.usegalaxy.eu`)
    2. Provide the **username** (usually the email address) and the **password** on the Galaxy server
    3. Connect

4. Add the files to the FTP server by dragging/dropping them or right clicking on them and uploading them

    The FTP transfer will start. We need to wait until they are done.

5. Open the Upload menu on the Galaxy server
6. Click on **Choose FTP file** on the bottom
7. Select files to import into the history
8. Click on **Start**

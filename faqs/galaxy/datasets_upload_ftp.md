---
title: Upload many files (>10) via FTP
area: data upload
box_type: tip
layout: faq
contributors: [bebatut, shiltemann]
---

Some Galaxies offer FTP upload for very large datasets.

**Note:** the *"Big Three"* Galaxies (Galaxy Main, Galaxy EU, and Galaxy Australia) no longer support FTP upload, due to the recent improvements
of the default web upload, which should now support large file uploads and almost all use cases. For situations where uploading via the web
interface is too tedious, the
[galaxy-upload commandline utility](https://github.com/galaxyproject/galaxy-upload) is also available as an alternative to FTP.

To upload files via FTP, please

1. Check that your Galaxy supports FTP upload and look up the FTP settings.

2. Make sure to have an FTP client installed

    There are many options. We can recommend [FileZilla](https://filezilla-project.org/), a free FTP client that is available on Windows, MacOS, and Linux.

3. Establish FTP connection to the Galaxy server
    1. Provide the Galaxy server's FTP server name (e.g. `ftp.mygalaxy.com`)
    2. Provide the **username** (usually the e-mail address) and the **password** on the Galaxy server
    3. Connect

4. Add the files to the FTP server by dragging/dropping them or right clicking on them and uploading them

    The FTP transfer will start. We need to wait until they are done.

5. Open the Upload menu on the Galaxy server
6. Click on **Choose FTP file** on the bottom
7. Select files to import into the history
8. Click on **Start**

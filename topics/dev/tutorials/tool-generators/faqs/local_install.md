---
title: Why do I need that big (~5GB!) complicated Docker thing - can I just install the ToolFactory into our local galaxy server from the toolshed?
box_type: question
layout: faq
contributors: [fubar2]
---

You can but it can’t really be very useful. The ToolFactory is a Galaxy tool, but it installs newly generated tools automatically into the local Galaxy server. This is not normally possible because a tool cannot escape Galaxy's job execution environment isolation. The ToolFactory needs to write to the normally forbidden server's configuration so the new tool appears in the tool menu and is installed in the TFtools directory which is a subdirectory of the Galaxy tools directory.
The Appliance is configured so the ToolFactory and the Planemo test tool use remote procedure calls (RPC using rpyc) to do what tools cannot normally do. The rpyc server runs in a separate container. Without it, tool installation and testing are difficult to do inside Galaxy tools.
Known good tools can be uploaded to a local toolshed from your private appliance for installation to that server of yours. Debugging tools on a production server is not secure SOP. You just never know what might break. *That’s why a desktop disposable appliance is a better choice.*


---
title: Why isn't my history updating?
layout: faq
area: debugging
box_type: tip
google_form_id: 1729851011
contributors:
- hexylena
---
Have you ever experienced that you would submit a job but your history wouldn't update? Maybe it doesn't scroll or the datasets stay permanently grey even when you know they should be complete, *until you refresh the webpage*?

One possible cause of this can be a difference in the clocks of your browser and the server. Check that your clocks match, and if not, reconfigure them! If you are following the Galaxy Admin Training, you will have setup `chrony`. Check that your `chrony` configuration is valid and requesting time from a local pool.


```
# chronyc -n sources
210 Number of sources = 1
MS Name/IP address         Stratum Poll Reach LastRx Last sample               
===============================================================================
^? 169.254.169.123               0   7     0     -     +0ns[   +0ns] +/-    0ns

```

This command should return some valid sources. THe above shows an example of a time source that isn't working, 0ns is not a realistic office and LastRx is empty. Instead it should look more like::

```
# chronyc -n sources
210 Number of sources = 5
MS Name/IP address         Stratum Poll Reach LastRx Last sample               
===============================================================================
^? 169.254.169.123               0   6     0     -     +0ns[   +0ns] +/-    0ns
^? 178.239.19.58                 0   6     0     -     +0ns[   +0ns] +/-    0ns
^? 194.104.0.153                 2   6     1     0   +138us[ +138us] +/-   30ms
^? 45.138.55.61                  1   6     1     1   -103us[ -103us] +/- 3158us
^? 178.239.19.57                 2   6     1     1   -301us[ -301us] +/- 3240us
```

Here we see a number of sources, with more plausible offsets and non-empty LastRx.

If your time was misconfigured, you might now see something like:

```
# chronyc -n tracking
Reference ID    : B950F724 (185.80.247.36)
Stratum         : 2
Ref time (UTC)  : Tue Oct 22 09:44:29 2024
System time     : 929.234680176 seconds slow of NTP time
```

as chrony slowly adjusts the system clock to match NTP time.

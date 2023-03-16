---
layout: tutorial_hands_on

title: "Deploying Wireguard for private mesh networking"
zenodo_link: ""
questions:
  - What is wireguard?
  - When is it useful?
  - Is it right for me?
objectives:
  - Setup a wireguard mesh across a few nodes
time_estimation: "60m"
key_points:
  - Wireguard is incredibly easy to deploy, and very secure.
contributions:
  authorship:
  - hexylena
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
  - type: "none"
    title: "Three or more VMs (they can be tiny, 1 CPU, <1GB RAM)"
subtopic: cloud
tags:
  - wireguard
  - networking
---

In this tutorial we will briefly cover what [Wireguard](https://www.wireguard.com/) is and how you can leverage it for your needs. This will not make you an expert on Wireguard but will give you the tools you need in order to setup a local Wireguard network.

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# What is Wireguard?

Wireguard is a VPN like OpenVPN or IPSec, but instead of the hub and spoke model of those where all traffic must go through a central node, Wireguard creates a mesh network where machines can all talk individually to each other. Wireguard also uses modern encryption only, ensuring your data stays safe.

## Is it right for me?

If you have machines that need to talk to each other privately, and you don't
have a better way to do it like a local network team, then yes, it's a great
solution to private, secure, fast networking.
It has [excellent performance](https://www.wireguard.com/performance/)
despite the encryption, and is built directly into the kernel.

By using wireguard, you can let services listen only on the wireguard interface, and thus only known and trusted machines can access those services.

# Setting up the infrastructure

> <hands-on-title>Configuration files</hands-on-title>
>
> 1. Create a `ansible.cfg` file (next to your playbook) to [configure settings](https://docs.ansible.com/ansible/2.9/reference_appendices/config.html) like the inventory file (and save ourselves some typing!), or the Python interpreter to use:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/ansible.cfg
>    @@ -0,0 +1,6 @@
>    +[defaults]
>    +interpreter_python = /usr/bin/python3
>    +inventory = hosts
>    +retry_files_enabled = false
>    +[ssh_connection]
>    +pipelining = true
>    {% endraw %}
>    ```
>    {: data-commit="Add ansible.cfg"}
>
>    > <tip-title>CentOS7</tip-title>
>    > As mentioned in the "Ubuntu or Debian, CentOS or RHEL?" comment above, if you are using CentOS7 do not set `interpreter_python` in `ansible.cfg` .
>    {: .tip}
>
>    Pipelining will make [Ansible run faster](https://docs.ansible.com/ansible/2.9/reference_appendices/config.html#ansible-pipelining) by significantly reducing the number of new SSH connections that must be opened.
>
> 2. Create the `hosts` inventory file if you have not done so yet, defining an `[wireguard]` group with every host you want to be part of the cluster. For each machine, also set a variable `wireguard_ip` with an address from `192.168.0.0/24`
>
>    > > <code-in-title>Bash</code-in-title>
>    > > ```bash
>    > > cat hosts
>    > > ```
>    > {: .code-in}
>    >
>    > > <code-out-title>Bash</code-out-title>
>    > >
>    > > Your hostname is probably different:
>    > >
>    > > {% raw %}
>    > > ```diff
>    > > --- /dev/null
>    > > +++ b/hosts
>    > > @@ -0,0 +1,5 @@
>    > > +[wireguard]
>    > > +1-wg.galaxy.training wireguard_ip=192.168.0.1
>    > > +2-wg.galaxy.training wireguard_ip=192.168.0.2
>    > > +3-wg.galaxy.training wireguard_ip=192.168.0.3
>    > > +4-wg.galaxy.training wireguard_ip=192.168.0.4
>    > > {% endraw %}
>    > > ```
>    > > {: data-commit="Add hosts"}
>    > {: .code-out}
>    {: .code-2col}
{: .hands_on}

Wireguard can use any of the [private network](https://en.wikipedia.org/wiki/Private_network) blocks, here we use `192.168.0.0/16` for familiarity, Tailscale uses `10.0.0.0/8`.

## Writing the playbook

First lets set up the playbook and install Wireguard, without configuring it.

> <hands-on-title>Installing Wireguard</hands-on-title>
>
> 1. Create and open `wg.yml` which will be our playbook. Add the following:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/wg.yml
>    @@ -0,0 +1,16 @@
>    +---
>    +- hosts: all
>    +  become: yes
>    +  vars:
>    +    wireguard_mask_bits: 24
>    +    wireguard_port: 51871
>    +  tasks:
>    +    - name: update packages
>    +      apt:
>    +        update_cache: yes
>    +        cache_valid_time: 3600
>    +
>    +    - name: Install wireguard
>    +      apt:
>    +        name: wireguard
>    +        state: present
>    {% endraw %}
>    ```
>    {: data-commit="Add initial playbook, install wireguard"}
>
> 2. Run the playbook:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook wg.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
{:.hands_on}

Now we can start configuring it. Wireguard relies on each node having a private/public keypair which is used for encrypting communications between the nodes. By default this uses modern elliptic cryptography which is quite simple and provably secure, wireguard only has ~4k LOC compared to 400k+ LOC for most VPN solutions. Let's generate those keys now:

> <hands-on-title>Generating a keypair</hands-on-title>
>
> 1. Edit `wg.yml` and add
>
>    {% raw %}
>    ```diff
>    --- a/wg.yml
>    +++ b/wg.yml
>    @@ -14,3 +14,8 @@
>           apt:
>             name: wireguard
>             state: present
>    +
>    +    - name: Generate Wireguard keypair
>    +      shell: wg genkey | tee /etc/wireguard/privatekey | wg pubkey | tee /etc/wireguard/publickey
>    +      args:
>    +        creates: /etc/wireguard/privatekey
>    {% endraw %}
>    ```
>    {: data-commit="Generate the keys"}
>
> 2. Run the playbook:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook wg.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Go check out what the keys look like!
{:.hands_on}

Ok, that's got wireguard setup, but it still isn't running. As of 2018 or so, systemd added built in support for wireguard and we'll use that to setup the network devices, and the network itself.

> <hands-on-title>Setting up the network</hands-on-title>
>
> 1. Edit `wg.yml` and add the following. We want to register the contents of these keys, as they'll be used to configure our networks later.
>
>    {% raw %}
>    ```diff
>    --- a/wg.yml
>    +++ b/wg.yml
>    @@ -19,3 +19,13 @@
>           shell: wg genkey | tee /etc/wireguard/privatekey | wg pubkey | tee /etc/wireguard/publickey
>           args:
>             creates: /etc/wireguard/privatekey
>    +
>    +    - name: register private key
>    +      shell: cat /etc/wireguard/privatekey
>    +      register: wireguard_private_key
>    +      changed_when: false
>    +
>    +    - name: register public key
>    +      shell: cat /etc/wireguard/publickey
>    +      register: wireguard_public_key
>    +      changed_when: false
>    {% endraw %}
>    ```
>    {: data-commit="Register the keys"}
>
> 2. Again editing `wg.yml`, we'll setup the network device:
>    {% raw %}
>    ```diff
>    --- a/wg.yml
>    +++ b/wg.yml
>    @@ -29,3 +29,28 @@
>           shell: cat /etc/wireguard/publickey
>           register: wireguard_public_key
>           changed_when: false
>    +
>    +    - name: Setup wg0 device
>    +      copy:
>    +        content: |
>    +          [NetDev]
>    +          Name=wg0
>    +          Kind=wireguard
>    +          Description=WireGuard tunnel wg0
>    +          [WireGuard]
>    +          ListenPort={{ wireguard_port }}
>    +          PrivateKey={{ wireguard_private_key.stdout }}
>    +          {% for peer in groups['all'] %}
>    +          {% if peer != inventory_hostname %}
>    +          [WireGuardPeer]
>    +          PublicKey={{ hostvars[peer].wireguard_public_key.stdout }}
>    +          AllowedIPs={{ hostvars[peer].wireguard_ip }}/32
>    +          Endpoint={{ hostvars[peer].inventory_hostname }}:{{ wireguard_port }}
>    +          PersistentKeepalive=25
>    +          {% endif %}
>    +          {% endfor %}
>    +        dest: /etc/systemd/network/99-wg0.netdev
>    +        owner: root
>    +        group: systemd-network
>    +        mode: 0640
>    +      notify: systemd network restart
>    {% endraw %}
>    ```
>    {: data-commit="Setup the netdev"}
>
>    Here we do a number of things:
>    - We configure a [NetDev](https://freedesktop.org/software/systemd/man/systemd.netdev.html), a virtual network device. Systemd supports many types but we'll use wireguard.
>    - Next we setup wireguard, specifying a `ListenPort` which is used for all wireguard communication, and specifies our PrivateKey which is used by the tunnel.
>    - For every peer (all other instances than ourselves), we setup a `WireGuardPeer` with that peer's public key, and which IPs are allowed to connect.
>
> 2. Again editing `wg.yml`, we'll setup the final bit, the network service.
>
>    {% raw %}
>    ```diff
>    --- a/wg.yml
>    +++ b/wg.yml
>    @@ -54,3 +54,16 @@
>             group: systemd-network
>             mode: 0640
>           notify: systemd network restart
>    +
>    +    - name: Setup wg0 network
>    +      copy:
>    +        content: |
>    +          [Match]
>    +          Name=wg0
>    +          [Network]
>    +          Address={{ wireguard_ip }}/{{ wireguard_mask_bits }}
>    +        dest: /etc/systemd/network/99-wg0.network
>    +        owner: root
>    +        group: systemd-network
>    +        mode: 0640
>    +      notify: systemd network restart
>    {% endraw %}
>    ```
>    {: data-commit="Setup the network"}
>
>    Here we just setup the network, declare the interface name, and the address for it.
>
> 2. And last let's restart the appropriate services.
>
>    {% raw %}
>    ```diff
>    --- a/wg.yml
>    +++ b/wg.yml
>    @@ -67,3 +67,10 @@
>             group: systemd-network
>             mode: 0640
>           notify: systemd network restart
>    +
>    +  handlers:
>    +    - name: systemd network restart
>    +      service:
>    +        name: systemd-networkd
>    +        state: restarted
>    +        enabled: yes
>    {% endraw %}
>    ```
>    {: data-commit="Restart the network"}
>
> 2. Run the playbook:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook wg.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Go check out the network! Try pinging each of the wireguard IPs to see if you can reach each of the machines.
{:.hands_on}

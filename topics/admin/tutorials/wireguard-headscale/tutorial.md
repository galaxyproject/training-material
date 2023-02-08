---
layout: tutorial_hands_on

title: "Deploying Tailscale/Headscale for private mesh networking"
zenodo_link: ""
questions:
  - What is Tailscale?
  - When is it useful?
  - Is it right for me?
objectives:
  - Setup a tailnet across a few nodes
time_estimation: "60m"
key_points:
  - Tailscale is a fantastic bit of software that Just Works™
  - We use headscale, an open source reimplementation of Tailscale's control server because it's easy to use in training
  - But if you can afford Tailscale, just use that.
  - There is a FOSS plan, go check it out!
contributions:
  authorship:
  - hexylena
  editing:
  - natefoo
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

[Tailscale](https://tailscale.com/) makes secure networking easy, it really is like magic. If you've used wireguard before, you know it takes a bit to setup and some configuration if you need to do anything fancy.


> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# What is Tailscale?

It's like Wireguard but easier and they built a **lot** of nice features on top of it. It's networking that "Just Works™" even more than Wireguard. If you prefer to use plain Wireguard without Headscale/Tailscale, or just want to get an understanding of the technology that Headscale/Tailscale build off of, [there is a tutorial for that]({% link topics/admin/tutorials/wireguard/tutorial.md %}) as well.

## Is it right for me?

if you have machines that need to talk to each other privately, and you don't
have a better way to do it like a local network team, then yes, it's a great
solution to private, secure, fast networking. if you need auditing, tailscale
will do that, rather than you having to build it out yourself.
it has [excellent performance](https://www.wireguard.com/performance/)
despite the encryption, and is built directly into the kernel.

By using wireguard, you can let services listen only on the wireguard
interface, and thus only known and trusted machines can access those services.

Tailscale makes wireguard setup even easier by removing the key management
step, which normally requires distributing keys to every machine. Instead that
step is handled centrally, and in the case of Tailscale enforceable with ACLs
and SSO and 2FA policies, however the networking remains meshed, and machines
connect directly to one another.

You can go one step further than trusting individual machines, with Tailscale,
as every device is tied to a single user, and you can ask Tailscale what is the
authenticated identity of the specific TCP connection, allowing automatically
logging in your users.

In the context of Galaxy, this can be useful for components like [Interactive 
Tools]({% link topics/admin/tutorials/interactive-tools/tutorial.md %}), which 
require a web proxy between Galaxy and the cluster node where the tool runs. If 
the cluster is not on the local network, Wireguard can be used to securely 
bridge the gap, and Headscale or Tailscale can greatly simplify that process.
# Setting up the infrastructure

{% include _includes/cyoa-choices.html option1="Tailscale" option2="Headscale" default="Tailscale" help="Here we provide a version of this tutorial that uses Headscale, if you want a more hands-on experience and feel like managing your own infrastructure." %}

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
> 2. Create the `hosts` inventory file if you have not done so yet.
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
>    > > <div class="Headscale" markdown="1">
>    > > Pick one to be the `head`, and the rest to be the `tail`. The head will act as the coordination server, the `tail` will be the nodes that should talk to each other.
>    > > {% raw %}
>    > > ```diff
>    > > --- /dev/null
>    > > +++ b/hosts
>    > > @@ -0,0 +1,5 @@
>    > > +[head]
>    > > +1-wg.galaxy.training
>    > > +[tail]
>    > > +2-wg.galaxy.training
>    > > +3-wg.galaxy.training
>    > > +4-wg.galaxy.training
>    > > {% endraw %}
>    > > ```
>    > > {: data-commit="Add hosts"}
>    > > </div>
>    > >
>    > > <div class="Tailscale" markdown="1">
>    > > Place all of the nodes in a `tail` group
>    > > ```diff
>    > > --- /dev/null
>    > > +++ b/hosts
>    > > @@ -0,0 +1,5 @@
>    > > +[tail]
>    > > +1-wg.galaxy.training
>    > > +2-wg.galaxy.training
>    > > +3-wg.galaxy.training
>    > > +4-wg.galaxy.training
>    > > ```
>    > > </div>
>    > {: .code-out}
>    {: .code-2col}
{: .hands_on}

<div class="Headscale" markdown="1">

> <details-title>Consider using Tailscale</details-title>
> Tailscale has implemented the coordination server and infrastructure with much better, more robust infrastructure that you won't have to be responsible for. If you can get your institution to pay for it, or grant money covering it, it's **probably worth it**. They add many new features often, and things like the mobile apps only work with their service. Tailscale is free for personal use (i.e. to test things out) and offers a free plan for open source projects that you may qualify for.
>
> For a training event obviously we want something free and quick to setup and destroy, so, we're using [Headscale](https://github.com/juanfont/headscale) since it's free and we're just going to destroy it immediately, and no one will accidentally get billed ;)
>
> Using Headscale will also teach you everything you need to know if you do choose to use Tailscale, which is simpler and has fewer components for you to manage yourself.
{: .details}

> <hands-on-title>Installing Headscale</hands-on-title>
>
> 1. Install the role
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles ckstevenson.headscale
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 2. Create and open `head.yml` which will be our playbook. Add the following:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/head.yml
>    @@ -0,0 +1,16 @@
>    +---
>    +- name: Headscale
>    +  hosts: head
>    +  become: true
>    +  vars:
>    +    headscale_user: 'headscale'
>    +    headscale_version: '0.15.0'
>    +    headscale_namespaces:
>    +    - galaxy
>    +  roles:
>    +    - ckstevenson.headscale
>    +  post_tasks:
>    +    - command: headscale --namespace galaxy preauthkeys create --reusable --expiration 1h
>    +      register: authkey
>    +    - debug:
>    +        msg: "{{ authkey.stdout.split('\n')[-1] }}"
>    {% endraw %}
>    ```
>    {: data-commit="Setup headscale"}
>
> 2. Run the playbook:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook headscale.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. This will return a code in the debug output. Save this code, you'll need it shortly
{:.hands_on}
</div>

Now we can setup the nodes

> <hands-on-title>Configure the nodes</hands-on-title>
>
> 1. Install the role
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles artis3n.tailscale
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 1. Edit `tail.yml` and add the following.
>
>    <div class="Headscale" markdown="1">
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/tail.yml
>    @@ -0,0 +1,17 @@
>    +---
>    +- name: Tailscale
>    +  hosts: tail
>    +  become: true
>    +  vars:
>    +    tailscale_args: "--advertise-exit-node --login-server http://{{ hostvars[groups['head'][0]].inventory_hostname }}:8080"
>    +  pre_tasks:
>    + - sysctl:
>    +     name: net.ipv4.ip_forward
>    +     value: '1'
>    +     state: present
>    + - sysctl:
>    +     name: net.ipv6.conf.all.forwarding
>    +     value: '1'
>    +     state: present
>    +  roles:
>    +    - artis3n.tailscale
>    {% endraw %}
>    ```
>    {: data-commit="Setup the clients"}
>    </div>
>
>    <div class="Tailscale" markdown="1">
>    ```diff
>    --- /dev/null
>    +++ b/tail.yml
>    @@ -0,0 +1,17 @@
>    +---
>    +- name: Tailscale
>    +  hosts: tail
>    +  become: true
>    +  vars:
>    +    tailscale_args: "--advertise-exit-node"
>    +  pre_tasks:
>    + - sysctl:
>    +     name: net.ipv4.ip_forward
>    +     value: '1'
>    +     state: present
>    + - sysctl:
>    +     name: net.ipv6.conf.all.forwarding
>    +     value: '1'
>    +     state: present
>    +  roles:
>    +    - artis3n.tailscale
>    ```
>    </div>
>
> 2. Run the playbook:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook tail.yml -e tailscale_authkey=YOUR_CODE
>    > ```
>    > {: data-cmd="true"}
>    >
>    > <div class="Headscale" markdown="1">
>    > Remember, you can find this code from the output of the first playbook
>    > </div>
>    > <div class="Tailscale" markdown="1">
>    > You can generate an authentication key under your [Tailscale account page](https://login.tailscale.com/admin/authkeys)
>    > </div>
>    {: .code-in}
>
> 3. Go check out your tailnet! Play around with the tailscale command and pinging other nodes with the suffix `.galaxy.example.com`{:.Headscale} `<username>.org.github.beta.tailscale.net`{:.Tailscale}
>
>    > Note that you'll need to enable [MagicDNS](https://login.tailscale.com/admin/dns) in your acount settings.
>    {:.Tailscale}
{:.hands_on}

> <tip-title>Exit Nodes</tip-title>
> We've configured `--advertise-exit-node`, which means you can direct ALL of your traffic to use one of your tailscale endpoints as an exit node, just run `tailscale up --exit-node=...`
>
> Note that:
> - If you're using headscale you need to manually enable that route (check the node list via `headscale nodes list` and then enable the specific route via `headscale nodes routes enable -i ...`), this is automatic in Tailscale
> - If you enable it, on a remote machine, it will immediately become unresponsive. Only do this on your local machine, e.g. a laptop connected to your tailnet.
{: .tip}

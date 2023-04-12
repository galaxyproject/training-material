---
layout: tutorial_hands_on

title: "Deploying a Beacon v1 in Galaxy"
zenodo_link: ""
questions:
    - What is a Beacon?
    - How do I deploy it?
    - Is v1 the same as v2?
objectives:
    - Deploy a Beacon
time_estimation: "30m"
key_points:
    - While deprecated, Beacon v1 is easy to deploy
    - It can also tick some boxes for grants!
contributions:
    authorship:
    - hexylena
    editing:
    - shiltemann
    funding:
    - CINECA-Project
subtopic: features
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
tags:
  - ga4gh
  - beacon
  - git-gat
---

This tutorial will guide you through setting up a [GA4GH Beacon](https://beacon-project.io/)!

> The Beacon Project is developed under a Global Alliance for Genomics and Health (GA4GH) Iniciative for the federated discovery of genomic data in biomedical research and clinical applications. 
{: .quote}

> <warning-title>Beacon v1</warning-title>
> This deploys an older Beacon v1 which was a simpler system.
> The Beacon v1 is more or less deprecated, with users being pushed to Beacon v2 which gives much richer answers, and offers better querying syntax.
{: .warning}

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Installing and Configuring

> <hands-on-title>Setting up a Beacon with Ansible</hands-on-title>
>
> 1. Setup the hosts
>
>    {% raw %}
>    ```diff
>    --- a/hosts
>    +++ b/hosts
>    @@ -6,3 +6,12 @@ galaxyservers
>     gat-0.au.training.galaxyproject.eu ansible_user=ubuntu
>     [monitoring]
>     gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    +
>    +[beacon]
>    +[beacon:children]
>    +beacon-import
>    +beacon-server
>    +[beacon-server]
>    +gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    +[beacon-import]
>    +gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    {% endraw %}
>    ```
>    {: data-commit="Add hosts"}
>
>    > <tip-title>What are these :children?</tip-title>
>    > Here we use some of the more advanced features of Ansible's Inventory system.
>    > We declare a host group called 'beacon' with no hosts of its own.
>    > 
>    > Then we declare that this beacon group has two children: beacon-import, and beacon server. We can then define host groups for those two entries with as many different hosts as we need. This makes it very easy to scale up configuration.
>    > 
>    > Here we will be using that feature to declare some 'beacon variables', which will be shared between the beacon-server and beacon-importer. Because they're children of 'beacon', they'll inherit any group variables defined for `group_vars/beacon.yml`.
>    {: .tip}
>
> 1. Setup the requirements
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -36,3 +36,7 @@
>       version: 2.1.3
>     - src: galaxyproject.proftpd
>       version: 0.3.1
>    +- name: paprikant.beacon
>    +  src: https://github.com/Paprikant/ansible-role-beacon
>    +- name: paprikant.beacon-importer
>    +  src: https://github.com/Paprikant/ansible-role-beacon_importer
>    {% endraw %}
>    ```
>    {: data-commit="Add requirements"}
>
> 2. Install the role with:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Create the vars file
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/group_vars/beacon.yml
>    @@ -0,0 +1,31 @@
>    +---
>    +postgres_data_dir: /data/beacon/postgresql/data
>    +postgres_init_dir: /data/beacon/postgresql/init
>    +bp_external_binding: 5050 # The default
>    +postgres_user: "{{ beacon_db_user }}"
>    +postgres_pass: "{{ beacon_db_password }}"
>    +postgres_external_binding: "{{ beacon_db_port }}"
>    +# Database Configuration
>    +beacon_db_user: beacon
>    +beacon_db_host: "{{ groups['beacon-server'][0] }}"
>    +beacon_db_password: "{{ vault_beacon_db_password }}"
>    +beacon_db_port: 9001
>    +#galaxy_api_key: This we will set in secrets.
>    +# Information about your beacon (consider filling this out.
>    +beacon_info_title: GA4GH Beacon
>    +beacon_info_beacon_id: your.galaxy.beacon
>    +beacon_info_description: Beacon service hosting datasets from all over the Galaxy
>    +beacon_info_url: https://{{ groups['beacon-server'][0] }}/beacon/
>    +beacon_info_service_group: galaxy-eu
>    +beacon_info_org_id: usegalaxy.aq
>    +beacon_info_org_name: Some Galaxy
>    +beacon_info_org_description: Galaxy community
>    +beacon_info_org_address: 123 Main Street, ZA
>    +beacon_info_org_welcome_url: https://galaxyproject.org/
>    +beacon_info_org_contact_url: https://galaxyproject.org/
>    +beacon_info_org_logo_url: https://galaxyproject.org/images/galaxy-logos/galaxy_project_logo_square.png
>    +beacon_info_org_info: More information about the organisation than just the description can go here.
>    +# Script Configuration
>    +galaxy_api_url: "https://{{ groups['galaxyservers'][0] }}"
>    +script_user: beacon
>    +script-dir: /home/beacon/script
>    {% endraw %}
>    ```
>    {: data-commit="Add relevant group variables"}
>
>    > <tip-title>groups?</tip-title>
>    > Here we again use some advanced features of Ansible's inventory system. Ansible knows the name of every hostname in the inventory. Now that we want to point the beacon configuration, either at the database which should be on `beacon-server`, or at Galaxy in `galaxyservers`, we ask the `groups` variable for what the inventory looks like. We use `[0]` to pull out the first hostname we find for both of those groups.
>    {: .tip}
>
> 3. Add the beacon-server playbook
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/beacon-server.yml
>    @@ -0,0 +1,9 @@
>    +---
>    +- name: Beacon Server
>    +  hosts: beacon-server
>    +  become: true
>    +  become_user: root
>    +  vars_files:
>    +    - group_vars/secret.yml
>    +  roles:
>    +    - paprikant.beacon
>    {% endraw %}
>    ```
>    {: data-commit="Add beacon server playbook"}
>
> 5. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook beacon-server.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in }
>
>    TODO: Check that it works
>
{: .hands_on}

## Setting up the Importer

Now that our beacon is running, we need to get data from Galaxy to the Beacon

> <hands-on-title>Setting up the Beacon Importer</hands-on-title>
>
> 1. Add the beacon-import playbook
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/beacon-import.yml
>    @@ -0,0 +1,9 @@
>    +---
>    +- name: Beacon Importer
>    +  hosts: beacon-import
>    +  become: true
>    +  become_user: root
>    +  vars_files:
>    +    - group_vars/secret.yml
>    +  roles:
>    +    - paprikant.beacon-importer
>    {% endraw %}
>    ```
>    {: data-commit="Add beacon importer playbook"}
>
> 1. Edit your `group_vars/secret.yml` and define some random passwords:
>
>    - The API key for your account, which must be an admin
>    
>    {% snippet faqs/galaxy/preferences_admin_api_key.md admin=true %}
>
>    > <code-in-title>Bash</code-in-title>
>    > ```
>    > ansible-vault edit group_vars/secret.yml
>    > ```
>    {: .code-in}
>
>    ```diff
>    --- a/group_vars/secret.yml
>    +++ b/group_vars/secret.yml
>    @@ -1,13 +1,15 @@
>    +galaxy_api_key: your-galaxy-api-key
>    +vault_beacon_db_password: some-super-secret-password
>    ```
>
>    <!-- Ignore this, just for the gat-automation. Vaults are ugly to work with :(
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/secret.yml
>    +++ b/group_vars/secret.yml
>    @@ -1,13 +1,17 @@
>     $ANSIBLE_VAULT;1.1;AES256
>    -62346261323266656232393034396134316636376533376139666437363535393562663838613938
>    -6336666266633563346337623265353935646361326337610a393834333233313461346439376438
>    -63383338346530656561636631666134373238366364363164313166346461383736613162653237
>    -3461363334323431370a656132303965653262386130353332623937376261396530393761353834
>    -38336565666437666436643163363831633331333766653266356163613138393734656465323634
>    -39366362383433366437353534663134313330316337393335383962613961386665633261616237
>    -35366635373063313631323939396164336330356361393464326636353037336461323531336434
>    -35613933303333623031353936393265636130363335376533393335663266313863376135383338
>    -36613464373231623938373434306266373234633036343636633963353361356631363533353066
>    -39323064336237646432323530313065303331326636353334343862373330313133326363363063
>    -38383564636161396435666164643334656435393533643163393434623434656238633631633939
>    -33353232666432376661
>    +62346462626130353731316233643332313965663232316430333231363531373038393233336439
>    +3234393435663964303162363933323762643039333732360a613338326139623964303265353330
>    +65303037653639373839316538626230623336313138346438393039633734383962343439653134
>    +3139663062666237300a633738383964393735393338313561306233643937313033303066363130
>    +32666266623766353838386434636130333263636563373739653739653562623834666135356234
>    +64643166366233656531656236396665643439626434353238326362326332626431323532396164
>    +35313762363333373463303936316233393033303239663238323739333133363362383935366562
>    +38376333663461326363633931663539313532376639373134313531663263386264636333623035
>    +32326637343630323265623062383962383963613231616230336238616562653039333964303262
>    +35643734626535376663366532633034616162396163626136613765666139613736303232336561
>    +36626337666362643339663935366232366632316662613466623235353934336635313063366139
>    +38376632656234313431313535626366393531636239626432343166633564353566356663343865
>    +39356133373134343235613332373331376135303636313232633664303539333962663535646561
>    +39643765356230313830633633396139333339613331343763323665326366383065316661636137
>    +38663138363633653434636433666665653632383964303433303538343630326337336438666530
>    +64313537623738616532
>    {% endraw %}
>    ```
>    {: data-commit="Add beacon/galaxy passwords to the vault"}
>
>    -->
>
> 1. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook beacon-importer.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in }
>
> 1. Add the nginx routes
>
>    {% raw %}
>    ```diff
>    --- a/templates/nginx/galaxy.j2
>    +++ b/templates/nginx/galaxy.j2
>    @@ -95,4 +95,14 @@ server {
>          location /reports/ {
>              proxy_pass http://unix:{{ galaxy_config.gravity.reports.bind }}:/;
>          }
>    +
>    +    location /beacon {
>    +        proxy_pass http://{{ groups['beacon-server'][0] }}:5050;
>    +        proxy_http_version 1.1;
>    +        proxy_set_header Upgrade $http_upgrade;
>    +        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
>    +        proxy_set_header X-Real-IP $remote_addr;
>    +        proxy_set_header Connection $connection_upgrade;
>    +        proxy_set_header Host $host;
>    +    }
>     }
>    {% endraw %}
>    ```
>    {: data-commit="Add beacon routes to nginx"}
>
> 1. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in }
{: .hands_on}

Congratulations, you've set up a Beacon v1 for Galaxy! Go check it out at `/beacon/` on your server.

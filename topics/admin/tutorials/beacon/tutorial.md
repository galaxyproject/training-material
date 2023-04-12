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
>    @@ -44,3 +44,7 @@
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
>    @@ -1,18 +1,22 @@
>     $ANSIBLE_VAULT;1.1;AES256
>    -39393239303639616131633130366134376431396464366131373430383435656261303633346638
>    -3965633961383235386230346561366637653561363961300a316537623964343132366163313038
>    -61663061656564393331663661643039386433326134636161636661613836396536636236336161
>    -3639636434333931380a646338383664646332343364393761313462346535623930663838313833
>    -31626164643766636564356333636364343164663562663733333261393039616535313439643264
>    -31333764393130326139306362333136656336663834316462333566313532373162353331393864
>    -37633064636364363335633331623063623639353536353263333661363166633636333833303862
>    -38343132313632336631356461323733616339646364626564323932373539343964373961313035
>    -31366532323539356532663537303763343035643066653062346238616233646439643061653537
>    -37636432343161323832336236326335626564393730663563353162303230306635393463643636
>    -33323165303061383732623066623837323037396563396263363230333634643739333262373932
>    -37633864383763616364316633366262666362643132623463346465373865373732346363353238
>    -64336230346636633637613265346630343231613161373861313932633638326431653534363933
>    -63383930313835653461633164386136343833646230656535356137626337373535373530376166
>    -38643930623339653739373132623731653536613135353236333235313137323135616636333239
>    -36663831653735396330366561346361623862626163386433626431353730323562626239373333
>    -37303031373237633934333236653361366563373332366166383662323733376335
>    +34633839393639313633633531313935666137366330303166643431323631343731353437383534
>    +3238386664303931383036326361306434336361613330650a326131353335373839643933333266
>    +35316266366231323830386633356631376132366433663133333233326134313662616636373432
>    +3164333962643366390a396531643764313439643666343062306433333665663162386335383366
>    +34313234393631393462636330343561313638343534613335333965373465663230393638316462
>    +34343639653863363565353239646138356664363737626638636635306132383765643335373930
>    +36623935663638613837633735666661343532316464383063323639343161376434383833353032
>    +33643332306337393530626138303364343037306335633630373433353635313264646565376166
>    +38336330616635323966313165373931313037356336383634623631323863336263356264623636
>    +61653062623539393433353461323031396539626263653735346331376638386437653639353063
>    +32303036336435663931393533373265326434646433643362316431303234633238383930303564
>    +64373631373861343738346562326232323661343362636432383661303866323932373565643734
>    +64616563313764643465376266626563616365643963373266353933356637316366326365366665
>    +31326531303065313835326662363330353466616665386565643865313331666530376336333430
>    +65303638346239343562656237343835656238626263393462653564333338616130323162333866
>    +33383662666432646131353233346531646533393630653031306364633165326364653337326262
>    +31343066323633366561613064316435306337393734386437313037376539316461353237316166
>    +36373263343038616631663866656531336339393063316365353162313235333463356339313131
>    +33386638636262636563633962356663666234326434663436383366313234346339623865373861
>    +62633732336438383564333234383763303763343737326338646263623138323265383261653634
>    +3336
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
>    @@ -104,4 +104,14 @@ server {
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

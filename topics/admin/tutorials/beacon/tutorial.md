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
>    +beacon_import
>    +beacon_server
>    +[beacon_server]
>    +gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    +[beacon_import]
>    +gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    {% endraw %}
>    ```
>    {: data-commit="Add hosts"}
>
>    > <tip-title>What are these :children?</tip-title>
>    > Here we use some of the more advanced features of Ansible's Inventory system.
>    > We declare a host group called 'beacon' with no hosts of its own.
>    > 
>    > Then we declare that this beacon group has two children: beacon_import, and beacon_server. We can then define host groups for those two entries with as many different hosts as we need. This makes it very easy to scale up configuration.
>    > 
>    > Here we will be using that feature to declare some 'beacon variables', which will be shared between the beacon_server and beacon_importer. Because they're children of 'beacon', they'll inherit any group variables defined for `group_vars/beacon.yml`.
>    {: .tip}
>
> 1. Setup the requirements
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -45,3 +45,7 @@
>       version: 2.1.5
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
>    +beacon_db_host: "{{ groups['beacon_server'][0] }}"
>    +beacon_db_password: "{{ vault_beacon_db_password }}"
>    +beacon_db_port: 9001
>    +#galaxy_api_key: This we will set in secrets.
>    +# Information about your beacon (consider filling this out.
>    +beacon_info_title: GA4GH Beacon
>    +beacon_info_beacon_id: your.galaxy.beacon
>    +beacon_info_description: Beacon service hosting datasets from all over the Galaxy
>    +beacon_info_url: https://{{ groups['beacon_server'][0] }}/beacon/
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
>    > Here we again use some advanced features of Ansible's inventory system. Ansible knows the name of every hostname in the inventory. Now that we want to point the beacon configuration, either at the database which should be on `beacon_server`, or at Galaxy in `galaxyservers`, we ask the `groups` variable for what the inventory looks like. We use `[0]` to pull out the first hostname we find for both of those groups.
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
>    +  hosts: beacon_server
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
>    +  hosts: beacon_import
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
>    @@ -1,22 +1,24 @@
>     $ANSIBLE_VAULT;1.1;AES256
>    -66653366633665383231303635643739396466653465626163633662623230643534346666303434
>    -3531356339646666623364383435306363373534373833370a333964306130653236373438373264
>    -34366433396135303064643932313135643064373034323865363939333565623238386161333465
>    -3932623162353034370a653663333363383938393936623663343639376430653366646662353830
>    -34383964326161383966306361633162653366353162633564623738353137333936313232363764
>    -38383962396637636136316366633936316365643038333030333932336365326335373737663834
>    -39363639303764343430336435313663313762353335366562656232646334356438343535303032
>    -35636131323163323036383163363964643237313137333131303737346662373233663562343031
>    -39396531306539623661376166313534303462623532323334393736316634616262313633626637
>    -37393066663239303564356232366466333334316565363631363230626665333039386133383933
>    -31656234363064646334623938383637393534313638643635396662643163346633346237353664
>    -34666161373764663166653737323736373261363565666633653831363833393264333666356633
>    -30373062333638366666316132623736386336383933326539663065653538386363323766313664
>    -63323737396264306439303366343834326164393763376663316366353766663131303966393039
>    -61343732313865613761623733356465336263633963376161663931653864383162393839643834
>    -34373239386461666462643162613536316166656135323332316163336635393035653164363138
>    -34336265383033663436306436663536383238366665653239313661666431666339323364646633
>    -64366361653034393962316562373835623631356339313062393933633634653735656636373031
>    -65376162633866623733623961303664663931636130373936326461623465363932373165356563
>    -36323838323865326464346337393164343765613333663830636234666663303535616532323936
>    -3230
>    +61363961643266653432653235396631323161383837383731633663653338326237346235646430
>    +3539643136643962316565333662323830386433633561390a363037333434626238366261323838
>    +66396362396135653836613933646637316361636464333166353034323261646630363065323933
>    +6132346534646439310a323364646531396366343837656164353262303833623935363137663535
>    +35623063333065383164643566353961373265326438323466633166626162333139306665633266
>    +31623130636639393361386262306237336132386563366435366631343038633462656661396638
>    +36356661363763656433303039623631613738626234393366636433323534353662636534306262
>    +36306461366635313861303539303739343465313161356433616634386632353338326338336539
>    +33326565396636363037313838383863346231303131396634666338303036393666353766326361
>    +36656461336562383861373839356564356163396337633939373337386661343638336362323663
>    +61613761343736613235316630366263343038313636393636323739313135643161653533346263
>    +63623633316432393738616635366236643962653530333134643064326631366365306166343235
>    +31656466373235323466653437306637643339656539383236393732383534333463396465663635
>    +38646561643733646231373230313066306237636630353763373035373832643763306566353636
>    +33346461623137316135633737373639633166313735333964656435646461323332646466343066
>    +30363232626566613561333965396361386166356630643061623564633766666631643738643663
>    +30303662313530396261633639653164326166353636633631386230646339616166353162313665
>    +36333162366136303532333636666339353831613663656463323666353230643834646434626666
>    +31366231303664396532613630616233653562336337653430396337356637633230303034656530
>    +33356662613636616634376163616638666164633565643135306166363435633034623935323438
>    +31383066633166663733626235313739383237376630643965333939313033386566623963316339
>    +31333264333939613866613136653137326339326537363466643966353938616361313934633566
>    +336130653132346566616361333861306131
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
>    +        proxy_pass http://{{ groups['beacon_server'][0] }}:5050;
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

{% snippet topics/admin/faqs/missed-something.md step=16 %}

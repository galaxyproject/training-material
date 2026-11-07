resource "openstack_compute_keypair_v2" "my-cloud-key" {
  name       = "my-key"
  public_key = "ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAACAQDV7gfNbNN5O8vH6/tM/iOFXKBP2YKRHXOmdfV8ogvu9BdVV0IPmDzk2EooVpThDE1VMv1hz3811tvBhHRJ6IgNhVIV/61w/+RazQD/AU27X8bX+Hb9EQ/bP4DW+6ySd/z5vdDLzpH5dbiMhzPEDkXVsylUT+hkQnas6cHspDhHmtKQ5MWOgDe3D/IEudTDJQe8hxxaU4TaZUmFzn7eYp9HvuK8qW0yCy4NWOxJJHA+G5wSCyLuKnaKo4AitUIzSKF1AB94oq7b96KONhPxgRptAk4OYIUTdNFbrI5HDaSNzHLnF5FbjQvG+Eu6m5nY5yvJMogE+jiuWeIXCZTCFljg287FUo0ohmbZpd802L6VXun14VumRC+rRgPrvBALo/CsyCsPIoBSTKhVElxKVOcRjmTLNfrUZM0GQxqJhIvah8BV+JTExkipPwkrKTdMAWIXvCoehxV+WMpBWqtEEzAzEoqJpaiec7HfriwsHTGESZWAPYEbFjzbHXQZtqBkbOvtokPMRmTWfWKxaplCMN6ddJeeY6faorD0w/e6lszWES1Q1ieajiPKDy37UvybKKvPTk4o3MzyzYOS4c8HQj+jnGeR5Q3ETuyz4psLyOfuBtIrfOeuxV42rFDmkYM3IrrRR+F9oklFG6Ig8DVfgQEzSG36NkgvpF4OdFvigYqXvw== cloud@vgcn"
}

resource "openstack_compute_instance_v2" "test" {
  name            = "central-manager"
  image_name      = "denbi-centos7-j10-2e08aa4bfa33-master"
  flavor_name     = "m1.tiny"
  key_pair        = "${openstack_compute_keypair_v2.my-cloud-key.name}"
  security_groups = ["default"]

  network {
    name = "public"
  }

  user_data = <<-EOF
    #cloud-config
    write_files:
    - content: |
        CONDOR_HOST = localhost
        ALLOW_WRITE = *
        ALLOW_READ = $(ALLOW_WRITE)
        ALLOW_NEGOTIATOR = $(ALLOW_WRITE)
        DAEMON_LIST = COLLECTOR, MASTER, NEGOTIATOR, SCHEDD
        FILESYSTEM_DOMAIN = terraform-training
        UID_DOMAIN = terraform-training
        TRUST_UID_DOMAIN = True
        SOFT_UID_DOMAIN = True
      owner: root:root
      path: /etc/condor/condor_config.local
      permissions: '0644'
    - content: |
        /data           /etc/auto.data          nfsvers=3
      owner: root:root
      path: /etc/auto.master.d/data.autofs
      permissions: '0644'
    - content: |
        share  -rw,hard,intr,nosuid,quota  ${openstack_compute_instance_v2.nfs.access_ip_v4}:/data/share
      owner: root:root
      path: /etc/auto.data
      permissions: '0644'
  EOF
}

resource "openstack_compute_instance_v2" "nfs" {
  name            = "nfs-server"
  image_name      = "denbi-centos7-j10-2e08aa4bfa33-master"
  flavor_name     = "m1.tiny"
  key_pair        = "${openstack_compute_keypair_v2.my-cloud-key.name}"
  security_groups = ["default"]

  network {
    name = "public"
  }

  user_data = <<-EOF
    #cloud-config
    write_files:
    - content: |
        /data/share *(rw,sync)
      owner: root:root
      path: /etc/exports
      permissions: '0644'
    runcmd:
     - [ mkdir, -p, /data/share ]
     - [ chown, "centos:centos", -R, /data/share ]
     - [ systemctl, enable, nfs-server ]
     - [ systemctl, start, nfs-server ]
     - [ exportfs, -avr ]
  EOF
}

resource "openstack_compute_instance_v2" "exec" {
  name            = "exec-${count.index}"
  image_name      = "denbi-centos7-j10-2e08aa4bfa33-master"
  flavor_name     = "m1.tiny"
  key_pair        = "${openstack_compute_keypair_v2.my-cloud-key.name}"
  security_groups = ["default"]
  count           = 2

  network {
    name = "public"
  }

  user_data = <<-EOF
    #cloud-config
    write_files:
    - content: |
        CONDOR_HOST = ${openstack_compute_instance_v2.test.access_ip_v4}
        ALLOW_WRITE = *
        ALLOW_READ = $(ALLOW_WRITE)
        ALLOW_ADMINISTRATOR = *
        ALLOW_NEGOTIATOR = $(ALLOW_ADMINISTRATOR)
        ALLOW_CONFIG = $(ALLOW_ADMINISTRATOR)
        ALLOW_DAEMON = $(ALLOW_ADMINISTRATOR)
        ALLOW_OWNER = $(ALLOW_ADMINISTRATOR)
        ALLOW_CLIENT = *
        DAEMON_LIST = MASTER, SCHEDD, STARTD
        FILESYSTEM_DOMAIN = terraform-training
        UID_DOMAIN = terraform-training
        TRUST_UID_DOMAIN = True
        SOFT_UID_DOMAIN = True
        # run with partitionable slots
        CLAIM_PARTITIONABLE_LEFTOVERS = True
        NUM_SLOTS = 1
        NUM_SLOTS_TYPE_1 = 1
        SLOT_TYPE_1 = 100%
        SLOT_TYPE_1_PARTITIONABLE = True
        ALLOW_PSLOT_PREEMPTION = False
        STARTD.PROPORTIONAL_SWAP_ASSIGNMENT = True
      owner: root:root
      path: /etc/condor/condor_config.local
      permissions: '0644'
    - content: |
        /data           /etc/auto.data          nfsvers=3
      owner: root:root
      path: /etc/auto.master.d/data.autofs
      permissions: '0644'
    - content: |
        share  -rw,hard,intr,nosuid,quota  ${openstack_compute_instance_v2.nfs.access_ip_v4}:/data/share
      owner: root:root
      path: /etc/auto.data
      permissions: '0644'
  EOF
}

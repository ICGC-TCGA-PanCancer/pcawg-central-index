use strict;
use Getopt::Long;

# TODO: 
# * need to include multiple run check
# * need to include setup of monitoring
# * need to include restart of failed workflows

my ($glob_target) = @ARGV;

my $test = 1;
my $verbose = 0;
my $glob_base = "";
my $glob_target = "target-*";
my $sensu_worker = "/glusterfs/netapp/homes1/BOCONNOR/gitroot/pancancer-sandbox/bionimbus_monitor/setup_sensu_worker.sh";
my $sensu_master = "/glusterfs/netapp/homes1/BOCONNOR/gitroot/pancancer-sandbox/bionimbus_monitor/setup_sensu_master.sh";

if (scalar(@ARGV) < 1 || scalar(@ARGV) > 6) {
 die "USAGE: perl $0 [--test] [--verbose] [--glob-base <path to directory that contains bindle dirs>] [--glob-target <target-*>]\n";
}

GetOptions(
  "test" => \$test,
  "verbose" => \$verbose,
  "glob-base=s" => \$glob_base,
  "glob-target=s" => \$glob_target, 
);

print "\n\n";
print "##############################\n";
print "# RUN DATE: ".`date`;
print "##############################\n";

my $glob_path = $glob_target;
if ($glob_base ne "") { $glob_path = "$glob_base/$glob_target"; }

foreach my $target (glob($glob_path)) {
  if (-d $target) {
    my $master_ip;
    #next if (defined($glob_target) && $glob_target ne '' && $glob_target ne $target );
    print "\n";
    print "##############################\n";
    print "# EXAMINING CLUSTER: $target #\n";
    print "##############################\n";
    foreach my $host ("$target/master", glob("$target/worker*")) {
      if (-d $host) {
        $host =~ /$target\/(\S+)/;
        my $hostname = $1;
        if ($hostname eq 'master') {
            $master_ip = `cd $host && vagrant ssh-config | grep HostName | awk '{print \$2}' 2> /dev/null`;
            chomp $master_ip;
        }
        print "\n##############################\n";
        print "#       HOST: $hostname      #\n";
        print "##############################\n\n";

        my $r = system("cd $host && vagrant ssh -c hostname 2> /dev/null");
        #print "  CMD STATUS: $r\n";
        # need to restart, ssh doesn't work!
        if ($r != 0) {
          print "  REBOOTING HOST\n";
          my $ip = `cd $host && vagrant ssh-config | grep HostName | awk '{print \$2}' 2> /dev/null`;
          chomp $ip;
          if ($ip =~ /\d+\.\d+\.\d+\.\d+/) {
            my $nova_id = `bash -l -c 'nova list | grep "$ip " | awk "{print \\\$2}"'`;
            chomp $nova_id;
            #print "  IP: $ip NOVA_ID $nova_id\n";
            # reboot here
            reboot_host($nova_id);
          }
        } else {
          # now check the hostname
          my $remote_name = `cd $host && vagrant ssh -c hostname 2> /dev/null`;
          #print "  REMOTE NAME: $remote_name\n";
          chomp $remote_name;
          $remote_name =~ /(\S+)/;
          $remote_name = $1;
          if ($remote_name ne $hostname) {
            print "  RESETTING HOST: REMOTE: $remote_name LOCAL: $hostname\n";
            reset_host($host, $hostname, $master_ip);
          } else {
            print "  HOST OK\n";
          }
        }
      }
    }
  }
}

sub reset_host {
  my ($dir, $host, $master_ip) = @_;
  #my $cmd = "cd $dir && vagrant ssh -c 'sudo hostname $host && sudo mount $master_ip:/home /home && sudo mount $master_ip:/mnt/home /mnt/home && sudo mount $master_ip:/mnt/datastore /mnt/datastore && sudo /etc/init.d/gridengine-exec restart; if [ -e /etc/init.d/gridengine-master ]; then sudo /etc/init.d/gridengine-master restart; fi;'";
  # FIXME: for now this will just disable hosts in SGE 
  my $cmd = "cd $dir && vagrant ssh -c 'sudo /etc/init.d/gridengine-exec stop; sudo hostname $host && sudo mount $master_ip:/home /home && sudo mount $master_ip:/mnt/home /mnt/home && sudo mount $master_ip:/mnt/datastore /mnt/datastore && sudo bash $sensu_worker' 2> /dev/null";
  if ($host eq 'master') {
    # FIXME: need to restart the master for sure
    $cmd = "cd $dir && vagrant ssh -c 'sudo /etc/init.d/gridengine-exec stop && sudo /etc/init.d/gridengine-master stop && sudo hostname $host && sudo /etc/init.d/gridengine-master start && sudo bash $sensu_master' 2> /dev/null";
  }
  if ($verbose) { print "    RESETTING CMD: $cmd\n"; }
  if (!$test) {
    my $r = system($cmd);
    if ($r) { print "PROBLEMS RESETTING HOST NAME\n"; }
  }
}

sub reboot_host {
  my ($nova_id) = @_;
  my $cmd = "bash -l -c 'nova reboot $nova_id'";
  if ($verbose) { print "    REBOOTING CMD: $cmd\n"; }
  if (!$test) {
    my $r = system($cmd);
    if ($r != 0) { print "PROBLEMS REBOOTING HOST\n"; }
  }
}


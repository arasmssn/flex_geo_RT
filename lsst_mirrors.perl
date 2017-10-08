#!/usr/bin/perl -w
# script to simulate optics for LSST
# motivated by the desire to use transformation directives in 
# the specifications for each optic

use Getopt::Long;
use Carp;
sub make_command ($$$);


my ($model,$band,$doparts,$filtermodel,$help);

($model,$band,$doparts,$filtermodel,$help)=("0904","I","BTC",undef,0);
($model,$band,$doparts,$filtermodel,$help)=("1008","I","BTC",undef,0);
undef($filtermodel);
# note - BTC is "baffles" - "telescope" - "camera" in any order

printf STDERR "ARGV : %s\n",join(' ',@ARGV);

usage() if (! GetOptions('model=s'  => \$model,
			 'band=s'   => \$band,
			 'filtermodel=s'  => \$filtermodel,
			 'doparts=s'=> \$doparts,
			 'help|?'   => \$help));
		  usage() if ($help);

{
    $band = lc $band;
    $doparts = lc $doparts;
    if (! defined($filtermodel)) {
	$filtermodel = $band;
    }
    $filtermodel = lc $filtermodel;

    printf STDERR "model $model band $band filtermodel $filtermodel doparts $doparts\n";

# check parameters to see if blessed values are being used.

# allowed models
    for my $mod ("1008","0904","0707","0608") {
	goto GOOD_MODEL if ($mod eq $model);
    }
    usage($model);
  GOOD_MODEL:
# allowed bands
    for my $bnd ("u","g","r","i","z","y2","y3","y4") {
	goto GOOD_BAND if ($bnd eq $band)
    }
    usage($band);
  GOOD_BAND:
    $z=0;

    setup_media();
    setup_models();

    $bafflecmds=[];
    $telescopecmds=[];
    $cameracmds=[];

    $cmdlist=$bafflecmds;

    # do the baffle & spider.

    @elements=();
    if ($doparts =~ /b/) {
	push(@elements,"CENTERCAP","CENTERPLATE","UPR_TEA_BFL","LWR_TEA_BFL",
	     "M2_BACK","M2_BAFFLE_BOTTOM","MRR_CVR_BFL","M1_APP_STOP");
    }

    foreach $element (@elements) {
	if (defined($$element{$model})) {
	    *model_element = $$element{$model};
	    if ($element eq "CENTERPLATE") {
		push(@{$cmdlist},
		     sprintf("baffle -T -t 0 0 %f -- ".
			     "-sn 4 -sw 800 -sl 3200 ".
			     "-ri %f -ro %f -rg 1",
			     $model_element{zo},
			     $model_element{ri},
			     $model_element{ro}),
		     sprintf("baffle -T -t 0 0 %f -d -r 0 0 45 -- ".
			     "-sn 4 -sw 1676.4 -sl 2845.4 ".
			     "-ri %f -ro %f -rg 1",
			     $model_element{zo},
			     $model_element{ri},
			     $model_element{ro}));
	    } else {
		if ($element eq "UPR_TEA_BFL") {
		    # include the spider/scissor etc.
		    my ($bfl_m,$bfl_z,$spl);
		    $bfl_m=1200/(4419.6-1600);
		    $bfl_z=$model_element{zo}-4419.6*$bfl_m;
		    $spl=($model_element{ri}+$model_element{ro});
		    # four separate spider components (actually 8 crosses)
		    # and finally the ring.
		    $include_spider=1;
		    if ($include_spider==1) {
			push(@{$cmdlist},
			     sprintf("baffle -T -t 0 +565.685 %f -- ".
				     "-sn 4 -sw 50 -sl %f -A1 %f -ri 0 -ro 0 -rg 1",
				     $bfl_z,$spl,$bfl_m),
			     sprintf("baffle -T -t 0 -565.685 %f -- ".
				     "-sn 4 -sw 50 -sl %f -A1 %f -ri 0 -ro 0 -rg 1",
				     $bfl_z,$spl,$bfl_m),
			     sprintf("baffle -T -t 0 +565.685 %f -- ".
				     "-sn 4 -sw 50 -sl %f -ri 0 -ro 0 -rg 1",
				     $model_element{zo},$spl),
			     sprintf("baffle -T -t 0 -565.685 %f -- ".
				     "-sn 4 -sw 50 -sl %f -ri 0 -ro 0 -rg 1",
				     $model_element{zo},$spl),
			     sprintf("baffle -T -t 0 0 %f -- ".
				     "-ri %f -ro %f -rg 1",
				     $model_element{zo},
				     $model_element{ri},
				     $model_element{ro}));
		    } else {
			# NO SPIDER!!
			push(@{$cmdlist},
			     sprintf("baffle -T -t 0 0 %f -- ".
				     "-ri %f -ro %f -rg 1",
				     $model_element{zo},
				     $model_element{ri},
				     $model_element{ro}));
		    }
		} else {
		    push(@{$cmdlist},
			 sprintf("baffle -T -t 0 0 %f -- ".
				 "-ri %f -ro %f -rg 1",
				 $model_element{zo},
				 $model_element{ri},
				 $model_element{ro}));
		}
	    }
	}
    }

    @elements=();

    push(@elements,"M1","M2CAM_HEXAPOD","M2","M3",
	 "CAMERA","L1","L2","FLT","L3","FP");

    foreach $element (@elements) {

	if ($element =~ /^M[123]$/) {
	    $cmdlist=$telescopecmds;
	} else {
	    $cmdlist=$cameracmds;
	}

	if (($element =~ /^M[123]/ or
	     $element =~ /^L[123]/ or
	     $element =~ /FLT/) and 
	    $element !~ /^M2CAM_HEXAPOD/) {

	    *model_element = $$element{$model};
	    *generic_switches = $$element{generic};
	    $extra="";
	    if ($element eq "M2") {
		$extra .= join(" ",@M2CAM_tran);
	    }
	    push(@{$cmdlist},join(' ',"asphere",
			       make_command(model_element,
					    generic_switches,
					    \$extra)));

	    if ($element =~ /FLT/ and 
		$band ne "i") {
		# correct current $z for the nonstandard filter thickness
		$z += (${$model_element{dz21}}{i} - 
		       ${$model_element{dz21}}{$band});
		printf STDERR "CORRECTING current z for nonstandard filter thickness %f\n",(${$model_element{dz21}}{i} - ${$model_element{dz21}}{$band});
	    }
#	    printf STDERR "command: %s\n",make_command(model_element);
	} else {

	    *model_element = $$element{$model};
	    *generic_switches = $$element{generic} if (defined($$element{generic}));

	    if ($element eq "CAMERA" or 
		$element eq "M2CAM_HEXAPOD") {
		# coordinate transform for general camera misalignments.
		# this is now kept in $$CAMERA{generic}{origin}.
		if ($element eq "M2CAM_HEXAPOD") {
		    # this needs to be handled as a "before & after" transform 
		    # for M2 and a "before only" transform before L1.
		    if (defined($generic_switches{origin}) and
			(defined($generic_switches{misalign}) and 
			 -e $generic_switches{misalign})) {
			my ($x0,$y0,$z0)=split(' ',$generic_switches{origin});
			push(@M2CAM_tran,"-T -t $x0 $y0 $z0 --");
			push(@M2CAM_tran,"`cat $generic_switches{misalign}`");
			push(@M2CAM_tran,"-T -t $x0 $y0 $z0 -I --");
#			$z -= $z0;
		    } else {
			@M2CAM_tran=();
		    }
		} else {
		    # first transform (non-inline) into the coordinate frame
		    # for the M2CAM hexapod. From there transform according to
		    # any CAM misalignments, if any.
		    # this requires a one-way coordinate transform as done 
		    # prior to interception by M2.

		    if ($#M2CAM_tran > -1) {
			printf STDERR "M2CAM_TRAN: %s\n",join(" ",@M2CAM_tran);
			foreach my $tran ( @M2CAM_tran ) {
			    push(@{$cmdlist},sprintf("tran_ray %s",$tran));
			}
		    }
		    # the camera misalignment, if any
		    if (defined($generic_switches{origin})) {
			my ($x0,$y0,$z0)=split(' ',$generic_switches{origin});

#			$z0 += 570.5747; # this would place origin near the FP

			printf STDERR "generic origin: %s\n",$generic_switches{origin};
			push(@{$cmdlist},sprintf("tran_ray -t $x0 $y0 $z0"));
			$z -= $z0;
		    }
		    # put in a coord transform to represent camera misalignment
		    # if $generic_switches{misalign}
		    if (defined($generic_switches{misalign}) and
			-e $generic_switches{misalign}) {
			push(@{$cmdlist},sprintf("tran_ray `cat %s`",
						 $generic_switches{misalign}));
		    }
		}
	    }

	    # simply propagate any shift to $z
	    $z += bandspec($model_element{dz});
	    if ($element eq "FP") {
		push(@{$cmdlist},sprintf("tran_ray -t 0 0 %f",$z));
		if (defined($generic_switches{misalign}) and
		    -e $generic_switches{misalign}) {
		    push(@{$cmdlist},sprintf("tran_ray `cat %s`",
					     $generic_switches{misalign}));
		}
		if (defined($generic_switches{solidbody_grav_distort}) and
		    -e $generic_switches{solidbody_grav_distort}) {
		    push(@{$cmdlist},sprintf("tran_ray `cat %s`",
					     $generic_switches{solidbody_grav_distort}));
		}
		push(@{$cmdlist},"plane_collapse -c 0 0 1 0");
	    }
	}
    }
    my @cmds=();
    push(@cmds,@{$bafflecmds})    if ($doparts =~ /b/);
    push(@cmds,@{$telescopecmds}) if ($doparts =~ /t/);
    push(@cmds,@{$cameracmds})    if ($doparts =~ /c/);
    $command=join(' | ',@cmds);
    undef(@cmds);
    printf STDERR "$command\n";
    open(F,"$command |") || die;
    while (<F>) {
	print;
    }
    open(FF,">cmds") || die;
    $command =~ s/\|/\n\|/g;
    printf FF "%s",$command;
    close(FF);
    exit;
}

sub make_command ($$$) {
    (*specs,*generic,*xtra)=@_;
#    printf STDERR "specs %s band %s\n",join(' ',%specs),$band;
    my @switch_list=();
    my ($lens)=(0);
    my (@arglist1,@arglist2);

    if (defined($specs{"R1"}) or defined($specs{"R2"}) or
	defined($specs{"z1"}) or defined($specs{"z2"}) or
	defined($specs{"k1"}) or defined($specs{"k2"})) {
#	printf STDERR "this looks like a lens..\n";
	$lens=1;
	push(@switch_list,"-L");
    }

    if ($lens==1) {
	@arglist1=("R1","k1","z1");
	@arglist2=("R2","k2","z2","ri","ro");
    } else {
	@arglist1=("R","k","z","ri","ro");
    }

    # the only ordering necessary is to name aspheric terms affecting surf.#1
    # directly after at least one surf.#1 parameter, and the same for surf.#2
    # in the case of a lens. There's no specification for sufrace #2 aspheric
    # terms except for in the ordering of those arguments.

    # use bandspec() as a dereferencer for values that are specified with a 
    # bandpass dependence.

    # use special parameters $specs{dz} (for mirrors and lenses)
    # and $specs{dz21} for lenses to account for the absolute axial 
    # positions of the apexes..
    
    $z=$z+bandspec($specs{dz});
    my ($transform_to_apex)=(1);

    push(@switch_list,$xtra);

    if ($transform_to_apex==1) {
	$specs{z}=0;
	$specs{z1}=0;
	my ($dx,$dy,$dz)=(0,0,$z);
	$transform_string=sprintf("-T -t $dx $dy $dz --");
	push(@switch_list,$transform_string);
    } else {
	$specs{z}=$z;
	$specs{z1}=$z;
    }

    foreach $key (@arglist1) {
	if (defined($specs{$key})) {
	    push(@switch_list,join(' ',"-".$key,bandspec($specs{$key})));
	}
    }
    
    if (($lens!=1) && defined($specs{asphs})) {
	*asph=$specs{asphs};
	foreach $key (sort keys %asph) {
	    push(@switch_list,"-".$key,bandspec($asph{$key}));
	}
    }
 
    if ($lens==1) {

	if ($transform_to_apex==1) {
	    $specs{z2}=bandspec($specs{dz21});
	} else {
	    $specs{z2}=$z+bandspec($specs{dz21});
	}

	$z += bandspec($specs{dz21});


	if (defined($specs{asphs1})) {
	    *asph=$specs{asphs1};
	    foreach $key (sort keys %asph) {
		push(@switch_list,"-".$key,bandspec($asph{$key}));
	    }
	}
	
	# add in 1st surface outside media specifications if listed (lenses)
	if (defined($generic{atm1})) {
	    $ref=bandspec($generic{atm1});
	    while (($key,$val) = each %{$ref}) {
		push(@switch_list,sprintf("-%s1 %s",$key,$val));
	    }
	}

	foreach $key (@arglist2) {
	    if (defined($specs{$key})) {
		push(@switch_list,join(' ',"-".$key,bandspec($specs{$key})));
	    }
	}

	if (defined($specs{asphs2})) {
	    *asph=$specs{asphs2};
	    foreach $key (sort keys %asph) {
		push(@switch_list,"-".$key,bandspec($asph{$key}));
	    }
	}
	# add in 2nd surface outside media specifications if listed (lenses)
	if (defined($generic{atm2})) {
	    $ref=bandspec($generic{atm2});
	    while (($key,$val) = each %{$ref}) {
		push(@switch_list,sprintf("-%s2 %s",$key,$val));
	    }
	}
    }

    # add in bandpass specifications if listed (filter)

    if (defined($generic{bandpass})) {
	if (defined($filtermodel)) {
	    $ref=${$generic{bandpass}}{$filtermodel};
	} else {
	    $ref=bandspec($generic{bandpass});
	}
	push(@switch_list,join(' ',"-I",@{$ref}));
    }

    # here add in solid body transformation directives if such files exist
    if (defined($generic{solidbody_misalign}) and 
	-e $generic{solidbody_misalign}) {
	push(@switch_list,sprintf("`cat %s`",$generic{solidbody_misalign}));
	
    }

    # here add in gravity distortions (solid body) transformation directives if such files exist
    if (defined($generic{solidbody_grav_distort}) and 
	-e $generic{solidbody_grav_distort}) {
	push(@switch_list,sprintf("`cat %s`",$generic{solidbody_grav_distort}));
	
    }

    # here add in zernike surface deformation directives if such files exist
    if (defined($generic{zernike_distortions}) and 
	-e $generic{zernike_distortions}) {
	push(@switch_list,sprintf("`cat %s`",$generic{zernike_distortions}));
	
    }
    
    # add in gravity induced surface deformation directives for each surface
    # if defined and if such files exist.

    if (defined($generic{zernike_grav_distort_s1}) and 
	-e $generic{zernike_grav_distort_s1}) {
	push(@switch_list,sprintf("`cat %s`",$generic{zernike_grav_distort_s1}));
    }

    if (defined($generic{zernike_grav_distort_s2}) and 
	-e $generic{zernike_grav_distort_s2}) {
	push(@switch_list,sprintf("`cat %s`",$generic{zernike_grav_distort_s2}));
    }

    # include distortion ntuple file if available
    if (defined($generic{fea_distortions_s1}) and 
	-e $generic{fea_distortions_s1}) {
	push(@switch_list,sprintf("-d1 %s",$generic{fea_distortions_s1}));
    }

    if (defined($generic{fea_distortions_s2}) and 
	-e $generic{fea_distortions_s2}) {
	push(@switch_list,sprintf("-d2 %s",$generic{fea_distortions_s2}));
    }

    # here add in FEA surface deformation directives if such files exist
    if (defined($generic{fea_distortions}) and 
	-e $generic{fea_distortions}) {
	push(@switch_list,sprintf("-d %s",$generic{fea_distortions}));
    }

    # here add in FEA surface deformation directives for each surface
    # if defined and if such files exist

    if (defined($generic{fea_distort_s1}) and 
	-e $generic{fea_distort_s1}) {
	push(@switch_list,sprintf("-d %s",$generic{fea_distort_s1}));
    }

    if (defined($generic{fea_distort_s2}) and 
	-e $generic{fea_distort_s2}) {
	push(@switch_list,sprintf("-d %s",$generic{fea_distort_s2}));
    }
    
    
    return(join(' ',@switch_list));
}

sub bandspec {
    local ($spec)=@_;
    return(${$spec}{$band}) if (defined(${$spec}{$band}));
    return($spec) if (defined($spec));
    return(0);
}


sub usage {
    print STDERR "Unknown option: @_\n" if ( @_ );
    print STDERR "usage: lsst_mirrors.perl [--model|-m MODEL] [--band|-b BAND] [--doparts|-d PARTS] [--help|-?]\n";
    print STDERR "       MODEL may be any of: (1008|0904|0707|0608|0511|0505|0407)\n";
    print STDERR "       BAND  may be any of: (U|G|R|I|Z|Y2|Y3|Y4)\n";
    print STDERR "       PARTS may be any combination of: {B,T,C} [baffles/telescope/camera]\n";
    exit;
}

sub setup_media {
    %atm_at_site=(p => 600.0/760.0,
		  t => 10,
		  h => 8);
    %atm_in_camera=(p => 600.0/760.0,
		    t => 15,
		    h => 0);
    %atm_in_cryostat=(p => 0,
		      t => 0,
		      h => 0);
    %atm_at_site=(p => 1,
		  t => 20,
		  h => 20);

    %atm_in_camera=%atm_at_site;
    %atm_in_cryostat=%atm_in_camera;

#    %lsst_bandpasses = (u => [3280*1.0250, 3980*1.0250],
#			g => [4100*1.0255, 5520*1.0255],
#			r => [5500*1.0255, 6940*1.0255],
#			i => [6940*1.0255, 8470*1.0255],
#			z => [8400*1.0255, 9500*1.0255],
#			y => [9600*1.0255,10280*1.0255]);

# new input based on collection-1777 (100819)

    %lsst_bandpasses = (u => [3210*1.0250, 3915*1.0250],
			g => [3997*1.0255, 5500*1.0255],
			r => [5490*1.0255, 6900*1.0255],
			i => [6880*1.0255, 8200*1.0255],
			z => [8185*1.0255, 9225*1.0255],
			y2=> [9790*1.0255,10245*1.0255],
			y3=> [9745*1.0255,11220*1.0255],
			y4=> [9225*1.0255,12000*1.0255]);
}

sub setup_models {
# setup the parameters for different models:

    %CENTERCAP  =("0904" => {ro => 1100,
                             ri => 0,
                             zo => -7418-2400});
    $CENTERCAP{"1008"}=$CENTERCAP{"0904"};

    %CENTERPLATE=("0904" => {ro => 0,
                             ri => 0,
                             zo => -7418-1200});
    $CENTERPLATE{"1008"}=$CENTERPLATE{"0904"};

    %UPR_TEA_BFL=("0904" => {ro => 5200.0,
                             ri => 4419.6,
                             zo => -7418});
    $UPR_TEA_BFL{"1008"}=$UPR_TEA_BFL{"0904"};

    %LWR_TEA_BFL=("0904" => {ro => 5125.0,
                             ri => 4350.0,
                             zo => -5377.2});
    $LWR_TEA_BFL{"1008"}=$LWR_TEA_BFL{"0904"};

    %MRR_CVR_BFL=("0904" => {ro => 4800,
                             ri => 4250,
                             zo => -3291.4});
    $MRR_CVR_BFL{"1008"}=$MRR_CVR_BFL{"0904"};

    %M2_BACK    =("0904" => {ro => 1735,
                             ri =>    0,
                             zo =>-6156.2007-300});

    %M2_BAFFLE_BOTTOM=("1008" => {ro => 2400,
				  ri =>    0,
				  zo =>-6156.2007-300+994});

    %M1_APP_STOP=("0904" => {ro => 4180,
                             ri => 4675,
                             zo => -439.4});
    $M1_APP_STOP{"1008"}=$M1_APP_STOP{"0904"};

    %M1=("generic" => {solidbody_misalign => "M1solidbody.dat",
		       solidbody_grav_distort => "M1solidbody_g.dat",
                       zernike_distortions=> "M1zern.dat",
		       fea_distortions => "M1_perturbation.tnt"},
         "0407" => {R     => -19200,
		    k     => -1.254809,
		    ro    => 4180,
		    ri    => 2573,
		    asphs => {A6 => 5.8773e-9}},
	 "0505" => {R     => -19400,
		    k     => -1.262,
		    ro    =>  4180,
		    ri    =>  2558,
		    asphs => {A6 => 7.940e-9}},
	 "0511" => {R     => -19835,
		    k     =>     -1.215,
		    ro    =>   4180,
		    ri    =>   2558,
                    dz    =>      0,
		    asphs => {A6 => 1.381e-9}},
	 "0608" => {R     => -19835,
		    k     =>     -1.218,
		    ro    =>   4180,
		    ri    =>   2558,
		    dz    =>      0,
		    asphs => {A6 => -0.406e-9,
			      A8 => 38.1e-12}}
	);

    $M1{"0707"}=$M1{"0608"};
    $M1{"0904"}=$M1{"0511"};
    $M1{"1008"}=$M1{"0904"};

# devise a coordinate transform that's applied "before & after" M2
# and then "before only" CAMSHIFT (and L1/L2/F/L3/FP)
# for the "before only" transform need extra bookkeeping because the "origin"
# for CAMSHIFT is written in absolute coordinates (WRT M1apex).

    # altered the origin from -7000 to +1075 so that planar regression
    # results are nearly insensitive to M2CAM rotations.
    %M2CAM_HEXAPOD=("generic" => {"misalign" => "M2CAMmisalign.dat",
                                  "origin"  => "0 0 1075"});
    # new origin based on fitting <fwhm> vs. origin for a fixed tilt of 0.01 deg
    %M2CAM_HEXAPOD=("generic" => {"misalign" => "M2CAMmisalign.dat",
                                  "origin"  => "0 0 -13480"});

    %M2=("generic" => {solidbody_misalign => "M2solidbody.dat",
		       solidbody_grav_distort => "M2solidbody_g.dat",
                       zernike_distortions=> "M2zern.dat",
		       fea_distortions => "M2_perturbation.tnt"},
         "0407" => {R     => -6032.68,
		    k     =>    -0.284514,
		    ro    =>  1600,
		    ri    =>   820,
		    asphs => {A6 => -1.2799e-5,
			      A8 => -9.1574e-7}},
	 "0505" => {R     => -6502.50,
		    k     => -0.43380,
		    ro    =>  1700,
		    ri    =>  900,
		    asphs =>  {A6 => -1.145e-5,
			       A8 => -5.689e-7}},
	 "0511" => {R     =>  -6788,
		    k     => -0.222,
		    ro    =>   1710,
		    ri    =>    900,
		    asphs => {A6 => -1.274e-5,
			      A8 => -9.68e-7}},
	 "0608" => {R     =>  -6789,
		    k     => -0.228,
		    ro    =>   1710,
		    ri    =>    900,
		    dz    =>  -6155.767,
		    asphs => {A6 => -1.380e-5,
			      A8 => -6.770e-7,
			      A10=> -3.230e-8}},
	 "0904" => {R     =>  -6788,
		    k     => -0.222,
		    ro    =>   1710,
		    ri    =>    900,
                    dz    =>  -6156.2007,
		    asphs => {A6 => -1.274e-5,
			      A8 => -9.68e-7}}

	 );

    $M2{"0707"}=$M2{"0608"};
    $M2{"1008"}=$M2{"0904"};

    %M3=("generic" => {solidbody_misalign => "M3solidbody.dat",
		       solidbody_grav_distort => "M3solidbody_g.dat",
                       zernike_distortions=> "M3zern.dat",
		       fea_distortions => "M3_perturbation.tnt"},
         "0407" => {R     => -8577.43,
		    k     =>     0.129492,
		    ro    =>  2573,
		    ri    =>  660,
		    asphs => {A6 => -2.2717e-7,
			      A8 => -3.6963e-9}},
	 "0505" => {R     => -8396.00,
		    k     =>  0.10120,
		    ro    =>  2478,
		    ri    =>  600,
		    asphs =>  {A6 => -3.353e-7,
			       A8 => -2.120e-9}},
	 "0511" => {R     => -8344.5,
		    k     =>     0.155,
		    ro    =>  2508,
		    ri    =>   527.35,
		    asphs => {A6 => -4.50e-7,
			      A8 => -8.15e-9}},
	 "0608" => {R     => -8344,
		    k     =>     0.155,
		    ro    =>  2508,
		    ri    =>   527.35,
		    dz    =>   6389.587,
		    asphs => {A6 => -4.47e-7,
			      A8 => -8.80e-9}},
	 "0904" => {R     => -8344.5,
		    k     =>     0.155,
		    ro    =>  2508,
		    ri    =>   550,
                    dz    => 6390,
		    asphs => {A6 => -4.50e-7,
			      A8 => -8.15e-9}}
	);

    $M3{"0707"}=$M3{"0608"};
    $M3{"1008"}=$M3{"0904"};

# camera origin is an arbitrary location, currently located about 572mm on the other side of the focal plane from M1/M3. Best if this were changed to the CG point for he camera. any misalignments of the camera will be measured from there - according to the tran_ray convention.
# UPDATE - origin of -32000 was chosen because this choice yielded focal surfaces that are (nearly) independent of rotations Rx & Ry.
    %CAMERA=("generic" => {misalign => "CAMmisalign.dat",
                           origin   => "0 0 -31850"},
             "0511" => {dz => {u => -3.505,
			       g => -1.7825,
			       r => -0.728,
			       i =>  0.00,
			       z =>  0.37,
			       y =>  0.5843,
			       y2=>  0.5843,
			       y3=>  0.5843,
			       y4=>  0.5843}},
             "0608" => {dz => {u => -3.383,
			       g => -1.833,
			       r => -0.658,
			       i =>  0.00,
			       z =>  0.374,
			       y =>  0.6848,
			       y2=>  0.6848,
			       y3=>  0.6848,
			       y4=>  0.6848}},
             "0904" => {dz => {u => -3.471,
			       g => -1.826,
			       r => -0.761,
			       i => -0.0993,
			       z =>  0.370,
			       y =>  0.5843,
			       y2=>  0.5843,
			       y3=>  0.5843,
			       y4=>  0.5843}},
             "1008" => {dz => {u => -3.5684,
			       g => -1.9166,
			       r => -0.7726,
			       i => -0.0821,
			       z =>  0.3225,
			       y =>  0.5719,
			       y2=>  0.5719,
			       y3=>  0.5719,
			       y4=>  0.5719}});
    ${$CAMERA{generic}}{origin}="0 0 -11225";

    $CAMERA{"0707"}=$CAMERA{"0608"};

    %L1=("generic" => {solidbody_misalign => "L1solidbody.dat",
		       solidbody_grav_distort  => "L1solidbody_g.dat",
		       fea_distortions_s1 => "L1_S1_perturbation_g.tnt",
		       fea_distortions_s2 => "L1_S2_perturbation_g.tnt",
		       zernike_grav_distort_s1 => "L1_S1_zern_g.dat",
		       zernike_grav_distort_s2 => "L1_S2_zern_g.dat",
                       zernike_distortions=> "L1zern.dat",
                       atm1               => \%atm_at_site,
                       atm2               => \%atm_in_camera},
         "0407" => {R1    => -2739.4,
		    R2    => -3803.2,
		    ro    =>   800,
		    ri    =>     0},
	 "0505" => {R1    => -2729.0,
		    R2    => -3722.0,
		    ro    =>   810,
		    ri    =>     0},
	 "0511" => {R1    => -2824,
		    R2    => -5021,
		    ro    =>   774,
		    ri    =>     0},
	 "0608" => {R1    => -2820,
		    R2    => -5017.5,
		    ro    =>   774,
		    ri    =>     0,
		    dz    => -3626.4,
		    dz21  =>   -82,
		    asphs1=> {A4 => -7.4e-6,
			      A6 => -3.0e-6}},
	 "0707" => {R1    => -2820,
		    R2    => -5017.5,
		    ro    =>   775,
		    ri    =>     0,
		    dz    => -3626.4,
		    dz21  =>   -82,
		    asphs1=> {A4 => -7.4e-6,
			      A6 => -3.0e-6}},
	 "0904" => {R1    => -2824,
		    R2    => -5021,
		    ro    =>   775,
		    ri    =>     0,
		    dz    => -3630.5,
		    dz21  =>   -82.23});

    $L1{"1008"}=$L1{"0904"};

    %L2=("generic" => {solidbody_misalign => "L2solidbody.dat",
		       solidbody_grav_distort => "L2solidbody_g.dat",
		       fea_distortions_s1 => "L2_S1_perturbation_g.tnt",
		       fea_distortions_s2 => "L2_S2_perturbation_g.tnt",
		       zernike_grav_distort_s1 => "L2_S1_zern_g.dat",
		       zernike_grav_distort_s2 => "L2_S2_zern_g.dat",
                       zernike_distortions=> "L2zern.dat",
		       atm1               => \%atm_in_camera,
		       atm2               => \%atm_in_camera},
         "0407" => {R1    => -5198.6,
		    R2    => -2058.5,
		    ro    =>   550,
		    ri    =>     0},
	 "0505" => {R1    => -6294.0,
		    R2    => -2038.0,
		    ro    =>   550,
		    ri    =>     0},
	 "0511" => {R1    =>     1e34,
		    R2    => -2529,
		    k1    =>     0,
		    k2    =>    -1.57,
		    ro    =>   551,
		    ri    =>     0,
		    asphs2=> {A6 => 1.656e-3}},
	 "0608" => {R1    =>     1e34,
		    R2    => -2528.4,
		    k1    =>     0,
		    k2    =>    -1.60,
		    ro    =>   551,
		    ri    =>     0,
		    dz    =>  -413.63,
		    dz21  =>   -30,
		    asphs2=> {A6 => 1.70e-3}},
	 "0707" => {R1    =>     1e34,
		    R2    => -2528.4,
		    k1    =>     0,
		    k2    =>    -1.60,
		    ro    =>   551,
		    ri    =>     0,
		    dz    =>  -413.63,
		    dz21  =>   -30,
		    asphs2=> {A6 => 1.70e-3}},
	 "0904" => {R1    =>     1e34,
		    R2    => -2529,
		    k1    =>     0,
		    k2    =>    -1.57,
		    ro    =>   551,
		    ri    =>     0,
		    dz    =>  -412.642,
		    dz21  =>   -30,
		    asphs2=> {A6 => 1.656e-3}}	 );

    $L2{"1008"}=$L2{"0904"};

    %FLT=("generic" => {solidbody_misalign => "Fsolidbody.dat",
			solidbody_grav_distort => "Fsolidbody_g.dat",
			fea_distortions_s1 => "F_S1_perturbation_g.tnt",
			fea_distortions_s2 => "F_S2_perturbation_g.tnt",
			zernike_grav_distort_s1 => "F_S1_zern_g.dat",
			zernike_grav_distort_s2 => "F_S2_zern_g.dat",
                        zernike_distortions=> "Fzern.dat",
                        atm1               => \%atm_in_camera,
                        atm2               => \%atm_in_camera,
                        bandpass           => \%lsst_bandpasses},
          "0407" => {R1   => -5630.32,
		     R2   => -5630.32,
		     ro   =>   390,
		     ri   =>     0},
	  "0505" => {R1   => -5986,
		     R2   => -5986,
		     ro   =>   385,
		     ri   =>     0},
	  "0511" => {R1    => -5624,
		     R2    => {u => -5507,
			       g => -5564,
			       r => -5597,
			       i => -5624,
			       z => -5624,
			       y => -5624,
			       y2=> -5624,
			       y3=> -5624,
			       y4=> -5624},
		     ro    => 375,
		     ri    =>   0},
	  "0608" => {R1    => -5620,
		     R2    => {u => -5514,
			       g => -5552,
			       r => -5585,
			       i => -5601,
			       z => -5611,
			       y => -5620,
			       y2=> -5620,
			       y3=> -5620,
			       y4=> -5620},
		     dz    =>  -360.37,
                     dz21  => {u => -26.3,
			       g => -21.5,
			       r => -17.8,
			       i => -15.7,
			       z => -14.5,
			       y => -13.5,
			       y2=> -13.5,
			       y3=> -13.5,
			       y4=> -13.5},
		     ro    => 375,
		     ri    =>   0},
	  "0707" => {R1    => -5620,
		     R2    => {u => -5514,
			       g => -5552,
			       r => -5585,
			       i => -5601,
			       z => -5611,
			       y => -5620,
			       y2=> -5620,
			       y3=> -5620,
			       y4=> -5620},
		     dz    =>  -360.37,
                     dz21  => {u => -26.3,
			       g => -21.5,
			       r => -17.8,
			       i => -15.7,
			       z => -14.5,
			       y => -13.5,
			       y2=> -13.5,
			       y3=> -13.5,
			       y4=> -13.5},
         	     ro    => 385,
		     ri    =>   0},
	  "0904" => {R1    => -5624,
		     R2    => {u => -5513, # check these: comparing recent
			       g => -5564, # gressler & chuck documents
			       r => -5594,
			       i => -5612,
			       z => -5624,
			       y => -5624,
			       y2=> -5624,
			       y3=> -5624,
			       y4=> -5624},
		     dz    =>  -357.58,
                     dz21  => {u => -26.2,
			       g => -21.14,
			       r => -17.8,
			       i => -15.7,
			       z => -14.2,
			       y => -13.5,
			       y2=> -13.5,
			       y3=> -13.5,
			       y4=> -13.5},
         	     ro    => 375,
		     ri    =>   0},
	  "1008" => {R1    => -5632,
		     R2    => {u => -5530, 
			       g => -5576, 
			       r => -5606,
			       i => -5623,
			       z => -5632,
			       y => -5640,
			       y2=> -5640,
			       y3=> -5640,
			       y4=> -5640},
		     dz    =>  -349.58,
                     dz21  => {u => -26.6,
			       g => -21.5,
			       r => -17.9,
			       i => -15.7,
			       z => -14.4,
			       y => -13.6,
			       y2=> -13.6,
			       y3=> -13.6,
			       y4=> -13.6},
         	     ro    => 375,
		     ri    =>   0}
	);

    # fudge all filters to be identical to $filtermodel band filter:
    if (0) {
	if (defined($filtermodel)) {
	    foreach my $bp ("u","g","r","i","z","y2","y3","y4") {
		next if ($bp eq $filtermodel);
		$FLT{"1008"}->{R2}->{$bp}=$FLT{"1008"}->{R2}->{$filtermodel};
		$FLT{"1008"}->{dz21}->{$bp}=$FLT{"1008"}->{dz21}->{$filtermodel};
	    }
	}
    }
    %L3=(generic => {solidbody_misalign => "L3solidbody.dat",
		     solidbody_grav_distort => "L3solidbody_g.dat",
		     fea_distortions_s1 => "L3_S1_perturbation_g.tnt",
		     fea_distortions_s2 => "L3_S2_perturbation_g.tnt",
		     zernike_grav_distort_s1 => "L3_S1_zern_g.dat",
		     zernike_grav_distort_s2 => "L3_S2_zern_g.dat",
                     zernike_distortions=> "L3zern.dat",
                     atm1               => \%atm_in_camera,
                     atm2               => \%atm_in_cryostat},
         "0407" => {R1    => -3625.5,
		    R2    => -17192.0,
		    ro    =>  365,
		    ri    =>  0	},
	 "0505" => {R1    => -2887.0,
		    R2    =>     1.0e34,
		    ro    =>   365,
		    ri    =>     0},
	 "0511" => {R1    => -3169,
		    R2    => 1.336e4,
		    k1    => -0.962,
		    ro    => 346,
		    ri    => 0},
	 "0608" => {R1    => -3142.4,
		    R2    => 1.386e4,
		    k1    => -0.937,
		    dz    =>    -45.3,
		    dz21  =>    -60,
		    ro    => 346,
		    ri    => 0},
	 "0707" => {R1    => -3142.4,
		    R2    => 1.386e4,
		    k1    => -0.937,
		    dz    =>   -45.3,
		    dz21  =>   -60,
		    ro    =>   365,
		    ri    =>     0},
	 "0904" => {R1    => -3169,
		    R2    => 1.336e4,
		    k1    => -0.962,
		    ro    => 361,
		    dz    =>    -45.3,
		    dz21  =>    -60,
		    ri    => 0},
	 "1008" => {R1    => -3169,
		    R2    => 1.336e4,
		    k1    => -0.962,
		    ro    => 361,
		    dz    =>    -53.3,
		    dz21  =>    -60,
		    ri    => 0}
	);

    %FP=("generic" => {misalign => "FPsolidbody.dat",
		       solidbody_grav_distort => "FPsolidbody_g.dat"},
         "0608" => {dz => -28.37},
         "0904" => {dz => -28.50});
    $FP{"0707"}=$FP{"0608"};
    $FP{"1008"}=$FP{"0904"};
}


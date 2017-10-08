#!/usr/bin/perl
# takes in the output of print_rays and outputs 
# a qdp file that is a summary of the distributions.

# printf "read serr 1 2\n";

$csi=0.8;
$wxmax=0.65;
$wymax=0.7;
$wxgap=0.02;
$wygap=0.03;
$gwxmin=0.1;
$gwxmax=0.9;
$gwymin=0.1;
$gwymax=0.9;
$convert_ccs=0;

while ($_=$ARGV[0],/^-/) {
    shift;
    /-c/ && (push(@comments,$ARGV[0]), shift, next);
    /-t/ && ($title=$ARGV[0], shift, next);
    /-r/ && ($convert_ccs=1, next);
}

printf "skip on\n";
printf "win 1\n";
printf "vi %f %f %f %f\n",$gwxmin,$gwymin,$wxmax,$wymax;
printf "ypl 1\n";
printf "ma 1 on 1\n";
printf "lab fi\n";
printf "!\n";
printf "win 2\n";
printf "vi %f %f %f %f\n",$gwxmin,$wymax+$wygap,$wxmax,$gwxmax;
printf "ypl 2\n";
printf "col 1 on 2\n";
printf "li on 2\n";
printf "lab nx off\n";
printf "!\n";
printf "win 3\n";
printf "vi %f %f %f %f\n",$wxmax+$wxgap,$gwymin,$gwxmax,$wymax;
printf "ypl 3\n";
printf "col 1 on 3\n";
printf "li on 3\n";
printf "lab ny off\n";
printf "!\n";
printf "mm 1 x 2\n";
printf "mm 1 y 3\n";
printf "!\n";

printf "win 1\n";
printf "lab x x\\dFP\\u [mm]\n";
printf "lab y y\\dFP\\u [mm]\n";
printf "fon rom\n";
printf "csi $csi\n";
printf "r y2 0\n";
printf "r x3 0\n";
printf "ti off\n";
printf "win all\nlab 50 vpos 0.5 0.95 csi 1.5 \"$title\"\n" if defined($title);

$labno=1;
$labypos=$gwymax;
foreach $com (@comments) {
    printf ("lab $labno col 1 vpos %f $labypos csi $csi jus lef \"$com\"\n",
	    $wxmax+$wxgap);
    $labno++;
    $labypos -= $wygap;
}

$xs=$ys=$n=0;
$xxs=$yys=0;

$xf="niceplot_xtmp";
$yf="niceplot_ytmp";
$xh="niceplot_xhtmp";
$yh="niceplot_yhtmp";

open(XF,">$xf");
open(YF,">$yf");

while (<>) {
    @_=split(' ');
#    print $_[4];
    if ($convert_ccs) {
	($x,$y)=@_[$[+8,$[+7]; 
    } else {
	# original
	($x,$y)=@_[$[+7,$[+8];
    }
    $xs+=$x;
    $ys+=$y;
    $xxs+=$x*$x;
    $yys+=$y*$y;
    $n++;
    printf XF "%s\n",$x;
    printf YF "%s\n",$y;
    printf "%s %s\n",$x,$y;
}

close(XF);
close(YF);

$sigmax=sqrt($xxs/$n-($xs/$n)*($xs/$n));
$sigmay=sqrt($yys/$n-($ys/$n)*($ys/$n));

$xmin=$xs/$n-5*$sigmax;
$xmax=$xs/$n+5*$sigmax;
$xn=300;
$ymin=$ys/$n-5*$sigmay;
$ymax=$ys/$n+5*$sigmay;
$yn=300;

printf "no no\n";
$histcom="histo 1 $xmin $xmax $xn < $xf > $xh";
print STDERR $histcom,"\n";
`$histcom`;
open(XF,"$xh");
while (<XF>) {
    if ($_ !~ /ead/) {
	($x,$xe,$y,$ye)=split(' ');
	printf ("%g %g\n",$x-$xe,$y);
	printf ("%g %g\n",$x+$xe,$y);
    }
}
close(XF);

printf "no no\n";
$histcom="histo 1 $ymin $ymax $yn < $yf > $yh";
print STDERR $histcom,"\n";
`$histcom`;
open(YF,"$yh");
while (<YF>) {
    if ($_ !~ /ead/) {
	($x,$xe,$y,$ye)=split(' ');
	printf ("%g %g\n",$y,$x-$xe);
	printf ("%g %g\n",$y,$x+$xe);
    }
}
close(YF);

$anacom="analyse_hist -c 3 $xh";
print STDERR $anacom,"\n";
print "win 2\n";
open(XF,"$anacom|");
while (<XF>) {
    if (/^50%/) {
	($wid,$cen)=($_ =~ /is\s(\S+)\scentered\sat\s(\S+)/);
	printf ("lab $labno col 2 pos %g 0 li 90 1 lst 2 \" \"\n",
		$cen-0.5*$wid);
	$labno++;
	printf ("lab $labno col 2 pos %g 0 li 90 1 lst 2 \" \"\n",
		$cen+0.5*$wid);
	$labno++;
    }
    if (/^90%/) {
	($wid,$cen)=($_ =~ /is\s(\S+)\scentered\sat\s(\S+)/);
	printf ("lab $labno col 4 pos %g 0 li 90 1 lst 2 \" \"\n",
		$cen-0.5*$wid);
	$labno++;
	printf ("lab $labno col 4 pos %g 0 li 90 1 lst 2 \" \"\n",
		$cen+0.5*$wid);
	$labno++;
    }
    if (/HEW/) {
	($fwhm,$hew)=($_ =~ /FWHM\s(\S+)\sHEW\s(\S+)\s/);
	printf ("lab fi (fwhm,hew)\\dx\\u=(%6.4f,%6.4f)\n",$fwhm,$hew);
    }
}
close(XF);

$anacom="analyse_hist -c 3 $yh";
print STDERR $anacom,"\n";
print "win 3\n";
open(YF,"$anacom|");
while (<YF>) {
    if (/^50%/) {
	($wid,$cen)=($_ =~ /is\s(\S+)\scentered\sat\s(\S+)/);
	printf ("lab $labno col 2 pos 0 %f li 0 1 lst 2 \" \"\n",
		$cen-0.5*$wid);
	$labno++;
	printf ("lab $labno col 2 pos 0 %f li 0 1 lst 2 \" \"\n",
		$cen+0.5*$wid);
	$labno++;
    }
    if (/^90%/) {
	($wid,$cen)=($_ =~ /is\s(\S+)\scentered\sat\s(\S+)/);
	printf ("lab $labno col 4 pos 0 %f li 0 1 lst 2 \" \"\n",
		$cen-0.5*$wid);
	$labno++;
	printf ("lab $labno col 4 pos 0 %f li 0 1 lst 2 \" \"\n",
		$cen+0.5*$wid);
	$labno++;
    }
    if (/HEW/) {
	($fwhm,$hew)=($_ =~ /FWHM\s(\S+)\sHEW\s(\S+)\s/);
	printf ("lab fi (fwhm,hew)\\dy\\u=(%6.4f,%6.4f)\n",$fwhm,$hew);
    }
}
close(YF);


package meshi.molecularElements.atoms;

import meshi.geometry.Coordinates;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomType;
import meshi.parameters.ResidueType;
import meshi.sequences.AtomAlignment;
import meshi.util.Rms;
import meshi.util.Utils;

public class AtomBuilder {
    //Based on pdb_from_torsions.pl (Below) by Dimaio lab.
    //https://faculty.washington.edu/dimaio/files/ 15.6.2019


    public static double C_N_len=1.33;
    public static double N_CA_len=1.45;
    public static double CA_C_len=1.52;
    public static double C_O_len=1.23;
    public static double CA_CB_len=1.52;
    public static double H_N_len=1.01;
    public static double HA_CA_len=1.09;
    public static double HB_CB_len=1.09;

    public static double C_N_CA_ang=58.3*Math.PI /180;
    public static double CA_C_O_ang=59.2*Math.PI /180;
    public static double N_CA_C_ang=68.8*Math.PI /180;
    public static double CA_C_N_ang=63.8*Math.PI /180;
    public static double N_CA_CB_ang=69.6*Math.PI /180;
    public static double H_N_CA_ang=60.8*Math.PI /180;
    public static double HA_CA_N_ang=71.9*Math.PI /180;
    public static double HA_CA_N_ang_gly=70.5*Math.PI /180;
    public static double HB_CB_CA_ang=70.5*Math.PI /180;
    public static double CB_OUT_OF_PLAIN = -122.8*Math.PI/180;
    /* Based on
    sub gen_atom_LA {
        my ($angle,$length,$X0,$X1) = @_;

        my $X10 = [$X0->[0]-$X1->[0],$X0->[1]-$X1->[1],$X0->[2]-$X1->[2]];
        my $X10p;
        if ($X10->[0] == 0 && $X10->[1] == 0) {
            $X10p=[0,1,0];
        } else {
            $X10p=[-$X10->[1],-$X10->[0],0];
        }
        my $W1 = cos( $angle/2 ); if ($angle<0) { $W1*=-1; }
        my $S1 = sqrt ( (1-$W1*$W1)/($X10p->[0]*$X10p->[0]+$X10p->[1]*$X10p->[1]+$X10p->[2]*$X10p->[2]) );
        my $R1 = quat2R( $X10p->[0]*$S1 , $X10p->[1]*$S1, $X10p->[2]*$S1, $W1 );
        my $RX10 = mapply($R1,$X10);
        normalize( $RX10 );

        my $Xnew = [$X0->[0]+$length*$RX10->[0],$X0->[1]+$length*$RX10->[1],$X0->[2]+$length*$RX10->[2]];
        return $Xnew;
    }*/
    public static void buildByAngle(double angle, double length, Atom atom0, Atom atom1, Atom target) {
        target.setXYZ(0,0,0);
        double[] X10 = {atom0.x() - atom1.x(), atom0.y() - atom1.y(), atom0.z() - atom1.z()};
        double[] X10p = new double[3];
        if ((X10[0] == 0) & (X10[1] == 0)) { X10p[0] = 0; X10p[1] = 1; X10p[2] = 0; }
        else { X10p[0] = -X10[1]; X10p[1] = -X10[0]; X10p[2] = 0; }
        double W1 = Math.cos(angle/2);
        if (angle <0)
            W1 *= -1;

        double S1 = Math.sqrt ( (1-W1*W1)/(X10p[0]*X10p[0]+X10p[1]*X10p[1]+X10p[2]*X10p[2]) );
        double[][] R1 = quat2R( X10p[0]*S1 , X10p[1]*S1, X10p[2]*S1, W1 );
        double[] RX10 = mapply(R1,X10);
        normalize( RX10 );
        target.setXYZ(atom0.x()+ length*RX10[0],atom0.y()+length*RX10[1],atom0.z()+length*RX10[2]);
    }

    /*
    Based on
    sub gen_atom_LAT {
        my ($torsion,$angle,$length,$X0,$X1,$X2) = @_;

        my $X10 = [$X0->[0]-$X1->[0],$X0->[1]-$X1->[1],$X0->[2]-$X1->[2]];
        my $X21 = [$X1->[0]-$X2->[0],$X1->[1]-$X2->[1],$X1->[2]-$X2->[2]];
        my $X10p = cross(  $X21, $X10 );
        normalize( $X10p );

        my $W1 = cos( $angle/2 ); if ($angle<0) { $W1*=-1; }
        my $S1 = sqrt ( (1-$W1*$W1)/($X10p->[0]*$X10p->[0]+$X10p->[1]*$X10p->[1]+$X10p->[2]*$X10p->[2]) );
        my $R1 = quat2R( $X10p->[0]*$S1 , $X10p->[1]*$S1, $X10p->[2]*$S1, $W1 );
        my $RX10 = mapply($R1,$X10);

        my $W2 = cos( $torsion/2 ); if ($torsion<0) { $W2*=-1; }
        my $S2 = sqrt ( (1-$W2*$W2)/($X10->[0]*$X10->[0]+$X10->[1]*$X10->[1]+$X10->[2]*$X10->[2]) );
        my $R2 = quat2R( $X10->[0]*$S2 , $X10->[1]*$S2, $X10->[2]*$S2, $W2 );
        my $RRX10 = mapply($R2,$RX10);
        normalize( $RRX10 );

        my $Xnew = [$X0->[0]+$length*$RRX10->[0],$X0->[1]+$length*$RRX10->[1],$X0->[2]+$length*$RRX10->[2]];

        return $Xnew;
    }
     */

    public static void buildByTorsion(double torsion, double angle, double length,Atom atom0,Atom atom1,Atom atom2, Atom target) {
        double[] X10 = {atom0.x()-atom1.x(),atom0.y()-atom1.y(),atom0.z()-atom1.z()};
        double[] X21 = {atom1.x()-atom2.x(),atom1.y()-atom2.y(),atom1.z()-atom2.z()};
        double[] X10p = cross(  X21, X10 );
        normalize( X10p );

        double W1 = Math.cos( angle/2 ); if (angle<0) { W1*=-1; }
        double S1 = Math.sqrt ( (1-W1*W1)/(X10p[0]*X10p[0]+X10p[1]*X10p[1]+X10p[2]*X10p[2]) );
        double[][] R1 = quat2R( X10p[0]*S1 , X10p[1]*S1, X10p[2]*S1, W1 );
        double[] RX10 = mapply(R1, X10);

        double W2 = Math.cos( torsion/2 ); if (torsion<0) { W2*=-1; }
        double S2 = Math.sqrt ( (1-W2*W2)/(X10[0]*X10[0]+X10[1]*X10[1]+X10[2]*X10[2]) );
        double[][] R2 = quat2R( X10[0]*S2 , X10[1]*S2, X10[2]*S2, W2 );
        double[] RRX10 = mapply(R2,RX10);
        normalize( RRX10 );
        target.setXYZ(atom0.x()+length*RRX10[0],atom0.y()+length*RRX10[1],atom0.z()+length*RRX10[2]);
    }



    /*
    Based on
    sub normalize {
        my $a = shift;
        my $b = sqrt($a->[0]*$a->[0] + $a->[1]*$a->[1] + $a->[2]*$a->[2]);
        if ($b > 1e-6) {
            $a->[0] /= $b; $a->[1] /= $b; $a->[2] /= $b;
        }
    }
    */
    public static void normalize(double[] a) {
        double b = Math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
        if (b > 1e-6) {
            a[0] /= b; a[1] /= b; a[2] /= b;
        }
    }

    /* Based on
        sub quat2R {
        my ($X,$Y,$Z,$W) = @_;
        my $xx = $X * $X; my $xy = $X * $Y; my $xz = $X * $Z;
        my $xw = $X * $W; my $yy = $Y * $Y; my $yz = $Y * $Z;
        my $yw = $Y * $W; my $zz = $Z * $Z; my $zw = $Z * $W;
        my $R = [ [ 1 - 2 * ( $yy+$zz ) ,     2 * ( $xy-$zw ) ,     2 * ( $xz+$yw ) ] ,
                  [     2 * ( $xy+$zw ) , 1 - 2 * ( $xx+$zz ) ,     2 * ( $yz-$xw ) ] ,
                  [     2 * ( $xz-$yw ) ,     2 * ( $yz+$xw ) , 1 - 2 * ( $xx+$yy ) ] ];
        return $R;
    }
     */
    public static double[][] quat2R( double X, double Y, double Z, double W) {
        double xx = X * X;
        double xy = X * Y;
        double xz = X * Z;
        double xw = X * W;

        double yy = Y * Y;
        double yz = Y * Z;
        double yw = Y * W;

        double zz = Z * Z;
        double zw = Z * W;
        double[][] R = new double[3][3];
        R[0][0] = 1 - 2 * ( yy+zz );  R[0][1] = 2 * ( xy-zw );     R[0][2] = 2 * ( xz+yw );
        R[1][0] = 2 * ( xy+zw );      R[1][1] = 1 - 2 * ( xx+zz ); R[1][2] = 2 * ( yz-xw );
        R[2][0] = 2 * ( xz-yw );      R[2][1] = 2 * ( yz+xw );     R[2][2] = 1 - 2 * ( xx+yy);
        return R;
    }

    /*
    Based on     sub mapply {
        my ($rotmat, $cart) = @_;
        my $out = [0, 0, 0];
        my ($i, $j);
        for ($i=0; $i < 3; ++$i) {
            for ($j=0; $j < 3; ++$j) {
                $out->[$i] += $rotmat->[$i][$j] * $cart->[$j];
            }
        }
        return $out;
    }
     */
    public static double[] mapply(double[][] rotmat, double[] cart) {
        double[] out = new double[3];
        for (int i=0; i < 3; i++) {
            for (int j=0; j < 3; j++) {
                out[i] += rotmat[i][j] * cart[j];
            }
        }
        return out;
    }

    /*
    Based on
        sub cross {
        my ($b,$c) = @_;
        my $a = [ $b->[1]*$c->[2] - $b->[2]*$c->[1] ,
                $b->[2]*$c->[0] - $b->[0]*$c->[2] ,
                $b->[0]*$c->[1] - $b->[1]*$c->[0] ];
        return $a;
    }
    */
    public static double[] cross(double[] b, double[] c){
        double[] a =  { b[1]*c[2] - b[2]*c[1] ,
                        b[2]*c[0] - b[0]*c[2] ,
                        b[0]*c[1] - b[1]*c[0] };
        return a;
    }


/*
    #!/usr/bin/perl
##
        ##
        ###############################################################################

    use strict;
    use constant PI    => 4 * atan2(1, 1);

if ($#ARGV < 0) {
        print STDERR "usage: $0 <topofile>\n";
        exit -1;
    }

    my $topofile = shift @ARGV;

    open (TOPO, $topofile) || die "Cannot open $topofile.";

    my @topolines;
while (<TOPO>) {
        my @f = split ' ', $_;
        push @topolines, \@f;
    }
    close(TOPO);

    my $C_N_len=1.33;
    my $N_CA_len=1.45;
    my $CA_C_len=1.52;
    my $C_O_len=1.23;
    my $CA_CB_len=1.52;
    my $H_N_len=1.01;
    my $HA_CA_len=1.09;
    my $HB_CB_len=1.09;

    my $C_N_CA_ang=58.3*PI/180;
    my $CA_C_O_ang=59.2*PI/180;
    my $N_CA_C_ang=68.8*PI/180;
    my $CA_C_N_ang=63.8*PI/180;
    my $N_CA_CB_ang=69.6*PI/180;
    my $H_N_CA_ang=60.8*PI/180;
    my $HA_CA_N_ang=71.9*PI/180;
    my $HA_CA_N_ang_gly=70.5*PI/180;
    my $HB_CB_CA_ang=70.5*PI/180;

    my ($AA,$N,$CA,$C,$O,$CB,$H1,$H2,$H3,$Oext,$HA1,$HA2,$HB1,$HB2,$HB3);
    my ($pN,$pCA,$pC);
    my $ctr=1;

# generate first residue
            $AA = "ALA";
if ($topolines[0][2] eq 'G') { $AA="GLY"; }
    $N=[0,0,0];
    $CA=gen_atom_L($N_CA_len,$N);
    $C=gen_atom_LA($N_CA_C_ang,$N_CA_len,$CA,$N);

    $H1=gen_atom_LAT(1/3*PI,$H_N_CA_ang,$H_N_len,$N,$CA,$C);
    $H2=gen_atom_LAT(PI,$H_N_CA_ang,$H_N_len,$N,$CA,$C);
    $H3=gen_atom_LAT(5/3*PI,$H_N_CA_ang,$H_N_len,$N,$CA,$C);

    $O=gen_atom_LAT(($topolines[0][1]*PI/180)-1*PI,$CA_C_O_ang,$C_O_len,$C,$CA,$N);
    printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " N  ", $AA, 1, $N->[0],$N->[1],$N->[2];
    printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " CA ", $AA, 1, $CA->[0],$CA->[1],$CA->[2];
    printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " C  ", $AA, 1, $C->[0],$C->[1],$C->[2];
    printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " O  ", $AA, 1, $O->[0],$O->[1],$O->[2];
if ($AA eq "ALA") {
        $CB=gen_atom_LAT((-122.8*PI/180),$N_CA_CB_ang,$CA_CB_len,$CA,$N,$C);
        printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " CB ", $AA, 1, $CB->[0],$CB->[1],$CB->[2];
        $HA1=gen_atom_LAT((-119*PI/180),$HA_CA_N_ang,$HA_CA_len,$CA,$N,$CB);
        printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " HA ", $AA, 1, $HA1->[0],$HA1->[1],$HA1->[2];
        $HB1=gen_atom_LAT(-1*PI,$HB_CB_CA_ang,$HB_CB_len,$CB,$CA,$N);
        $HB2=gen_atom_LAT(-1*PI/3,$HB_CB_CA_ang,$HB_CB_len,$CB,$CA,$N);
        $HB3=gen_atom_LAT(PI/3,$HB_CB_CA_ang,$HB_CB_len,$CB,$CA,$N);
        printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "1HB ", $AA, 1, $HB1->[0],$HB1->[1],$HB1->[2];
        printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "2HB ", $AA, 1, $HB2->[0],$HB2->[1],$HB2->[2];
        printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "3HB ", $AA, 1, $HB3->[0],$HB3->[1],$HB3->[2];
    } else {
        $HA1=gen_atom_LAT((121.4*PI/180),$HA_CA_N_ang_gly,$HA_CA_len,$CA,$N,$C);
        $HA2=gen_atom_LAT((117.2*PI/180),$HA_CA_N_ang_gly,$HA_CA_len,$CA,$N,$HA1);
        printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "1HA ", $AA, 1, $HA1->[0],$HA1->[1],$HA1->[2];
        printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "2HA ", $AA, 1, $HA1->[0],$HA1->[1],$HA1->[2];
    }
    printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "1H  ", $AA, 1, $H1->[0],$H1->[1],$H1->[2];
    printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "2H  ", $AA, 1, $H2->[0],$H2->[1],$H2->[2];
    printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "3H  ", $AA, 1, $H3->[0],$H3->[1],$H3->[2];

# additional
    foreach my $i (1..$#topolines) {
        $AA = "ALA";
        if ($topolines[$i][2] eq 'G') { $AA="GLY"; }

        $pN=$N;
        $pCA=$CA;
        $pC=$C;

        $N =gen_atom_LAT($topolines[$i-1][1]*PI/180,$CA_C_N_ang,$C_N_len,$pC,$pCA,$pN);
        $CA=gen_atom_LAT(PI,$CA_C_N_ang,$N_CA_len,$N,$pC,$pCA);
        $C =gen_atom_LAT($topolines[$i][0]*PI/180,$N_CA_C_ang,$CA_C_len,$CA,$N,$pC);
        $O =gen_atom_LAT($topolines[$i][1]*PI/180-1*PI,$CA_C_O_ang,$C_O_len,$C,$CA,$N);
        $H1=gen_atom_LAT(PI,$H_N_CA_ang,$H_N_len,$N,$CA,$pC);

        printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " N  ", $AA, $i+1, $N->[0],$N->[1],$N->[2];
        printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " CA ", $AA, $i+1, $CA->[0],$CA->[1],$CA->[2];
        printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " C  ", $AA, $i+1, $C->[0],$C->[1],$C->[2];
        printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " O  ", $AA, $i+1, $O->[0],$O->[1],$O->[2];
        if ($AA eq "ALA") {
            $CB=gen_atom_LAT((-122.8*PI/180),$N_CA_CB_ang,$CA_CB_len,$CA,$N,$C);
            printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " CB ", $AA, $i+1, $CB->[0],$CB->[1],$CB->[2];
            $HA1=gen_atom_LAT((-119*PI/180),$HA_CA_N_ang,$HA_CA_len,$CA,$N,$CB);
            printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " HA ", $AA, $i+1, $HA1->[0],$HA1->[1],$HA1->[2];
            $HB1=gen_atom_LAT(-1*PI,$HB_CB_CA_ang,$HB_CB_len,$CB,$CA,$N);
            $HB2=gen_atom_LAT(-1*PI/3,$HB_CB_CA_ang,$HB_CB_len,$CB,$CA,$N);
            $HB3=gen_atom_LAT(PI/3,$HB_CB_CA_ang,$HB_CB_len,$CB,$CA,$N);
            printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "1HB ", $AA, $i+1, $HB1->[0],$HB1->[1],$HB1->[2];
            printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "2HB ", $AA, $i+1, $HB2->[0],$HB2->[1],$HB2->[2];
            printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "3HB ", $AA, $i+1, $HB3->[0],$HB3->[1],$HB3->[2];
        } else {
            $HA1=gen_atom_LAT((121.4*PI/180),$HA_CA_N_ang_gly,$HA_CA_len,$CA,$N,$C);
            $HA2=gen_atom_LAT((117.2*PI/180),$HA_CA_N_ang_gly,$HA_CA_len,$CA,$N,$HA1);
            printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "1HA ", $AA, 1, $HA1->[0],$HA1->[1],$HA1->[2];
            printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , "2HA ", $AA, 1, $HA1->[0],$HA1->[1],$HA1->[2];
        }
        if ($i == $#topolines) {
            $Oext=gen_atom_LAT($topolines[$i][1]*PI/180,$CA_C_O_ang,$C_O_len,$C,$CA,$N);
            printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " Oxt", $AA, $i+1, $Oext->[0],$Oext->[1],$Oext->[2];
        }
        printf "ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", $ctr++ , " H  ", $AA, $i+1, $H1->[0],$H1->[1],$H1->[2];
    }

    sub gen_atom_LAT {
        my ($torsion,$angle,$length,$X0,$X1,$X2) = @_;

        my $X10 = [$X0->[0]-$X1->[0],$X0->[1]-$X1->[1],$X0->[2]-$X1->[2]];
        my $X21 = [$X1->[0]-$X2->[0],$X1->[1]-$X2->[1],$X1->[2]-$X2->[2]];
        my $X10p = cross(  $X21, $X10 );
        normalize( $X10p );

        my $W1 = cos( $angle/2 ); if ($angle<0) { $W1*=-1; }
        my $S1 = sqrt ( (1-$W1*$W1)/($X10p->[0]*$X10p->[0]+$X10p->[1]*$X10p->[1]+$X10p->[2]*$X10p->[2]) );
        my $R1 = quat2R( $X10p->[0]*$S1 , $X10p->[1]*$S1, $X10p->[2]*$S1, $W1 );
        my $RX10 = mapply($R1,$X10);

        my $W2 = cos( $torsion/2 ); if ($torsion<0) { $W2*=-1; }
        my $S2 = sqrt ( (1-$W2*$W2)/($X10->[0]*$X10->[0]+$X10->[1]*$X10->[1]+$X10->[2]*$X10->[2]) );
        my $R2 = quat2R( $X10->[0]*$S2 , $X10->[1]*$S2, $X10->[2]*$S2, $W2 );
        my $RRX10 = mapply($R2,$RX10);
        normalize( $RRX10 );

        my $Xnew = [$X0->[0]+$length*$RRX10->[0],$X0->[1]+$length*$RRX10->[1],$X0->[2]+$length*$RRX10->[2]];

        return $Xnew;
    }

    sub gen_atom_LA {
        my ($angle,$length,$X0,$X1) = @_;

        my $X10 = [$X0->[0]-$X1->[0],$X0->[1]-$X1->[1],$X0->[2]-$X1->[2]];
        my $X10p;
        if ($X10->[0] == 0 && $X10->[1] == 0) {
            $X10p=[0,1,0];
        } else {
            $X10p=[-$X10->[1],-$X10->[0],0];
        }
        my $W1 = cos( $angle/2 ); if ($angle<0) { $W1*=-1; }
        my $S1 = sqrt ( (1-$W1*$W1)/($X10p->[0]*$X10p->[0]+$X10p->[1]*$X10p->[1]+$X10p->[2]*$X10p->[2]) );
        my $R1 = quat2R( $X10p->[0]*$S1 , $X10p->[1]*$S1, $X10p->[2]*$S1, $W1 );
        my $RX10 = mapply($R1,$X10);
        normalize( $RX10 );

        my $Xnew = [$X0->[0]+$length*$RX10->[0],$X0->[1]+$length*$RX10->[1],$X0->[2]+$length*$RX10->[2]];
        return $Xnew;
    }

    sub gen_atom_L {
        my ($length,$X0) = @_;

        my $Xnew = [$X0->[0],$X0->[1],$X0->[2]+$length];

        return $Xnew;
    }

    sub cross {
        my ($b,$c) = @_;
        my $a = [ $b->[1]*$c->[2] - $b->[2]*$c->[1] ,
                $b->[2]*$c->[0] - $b->[0]*$c->[2] ,
                $b->[0]*$c->[1] - $b->[1]*$c->[0] ];
        return $a;
    }

    sub distance {
        my ($x1, $x2) = @_;
        my $i;
        my $out = sqrt (
                ($x1->[0]-$x2->[0])*($x1->[0]-$x2->[0]) +
                        ($x1->[1]-$x2->[1])*($x1->[1]-$x2->[1]) +
                        ($x1->[2]-$x2->[2])*($x1->[2]-$x2->[2])
        );
        return $out;
    }

    sub quatnorm {
        my ($X,$Y,$Z,$W) = @_;
        my $S = sqrt( $X*$X+$Y*$Y+$Z*$Z+$W*$W );
        return [ $X/$S , $Y/$S , $Z/$S , $W/$S ];
    }

    sub quat2R {
        my ($X,$Y,$Z,$W) = @_;
        my $xx = $X * $X; my $xy = $X * $Y; my $xz = $X * $Z;
        my $xw = $X * $W; my $yy = $Y * $Y; my $yz = $Y * $Z;
        my $yw = $Y * $W; my $zz = $Z * $Z; my $zw = $Z * $W;
        my $R = [ [ 1 - 2 * ( $yy+$zz ) ,     2 * ( $xy-$zw ) ,     2 * ( $xz+$yw ) ] ,
                  [     2 * ( $xy+$zw ) , 1 - 2 * ( $xx+$zz ) ,     2 * ( $yz-$xw ) ] ,
                  [     2 * ( $xz-$yw ) ,     2 * ( $yz+$xw ) , 1 - 2 * ( $xx+$yy ) ] ];
        return $R;
    }

    sub normalize {
        my $a = shift;
        my $b = sqrt($a->[0]*$a->[0] + $a->[1]*$a->[1] + $a->[2]*$a->[2]);
        if ($b > 1e-6) {
            $a->[0] /= $b; $a->[1] /= $b; $a->[2] /= $b;
        }
    }
    sub mapply {
        my ($rotmat, $cart) = @_;
        my $out = [0, 0, 0];
        my ($i, $j);
        for ($i=0; $i < 3; ++$i) {
            for ($j=0; $j < 3; ++$j) {
                $out->[$i] += $rotmat->[$i][$j] * $cart->[$j];
            }
        }
        return $out;
    }
    */
    //-------------------------------------- Add Side Chains --------------------------------------
    public static void addSideChains(Protein protein, String templateFile) {
        Protein template = new Protein(new AtomList(templateFile), ResidueExtendedAtoms.creator);
        for (Residue residue : protein.residues()) {
            if ((!residue.dummy()) & (residue.type != ResidueType.GLY) & (residue.type != ResidueType.ALA)) {
                Residue templateResidue = getTemplateResidue(template, residue.type());
                AtomList residueBBatoms = getBBAtoms(residue);
                AtomList templateBBatoms = getBBAtoms(templateResidue);
                AtomAlignment atomAlignment  = new AtomAlignment(residueBBatoms, templateBBatoms);
                Rms rms = new Rms(atomAlignment);
                Coordinates residueCOM = rms.getCenterOfMass0();
                Coordinates templateCOM = rms.getCenterOfMass1();
                for (Atom atom : templateResidue.getAtoms()) {
                    if ((atom != null) & (!atom.nowhere()))
                        atom.setXYZ(atom.x() - templateCOM.x(), atom.y() - templateCOM.y(), atom.z() - templateCOM.z());
                }
                for (Atom atom : residue.getAtoms()) {
                    if ((atom != null) & (!atom.nowhere()))
                        atom.setXYZ(atom.x() - residueCOM.x(), atom.y() - residueCOM.y(), atom.z() - residueCOM.z());
                }
                double[][] rotateMatrix = rms.getMatrix();
                for (Atom atom : templateResidue.getAtoms()) {
                    if (! atom.nowhere()) {
                        atom.setXYZ(
                                atom.x() * rotateMatrix[0][0] + atom.y() * rotateMatrix[0][1] + atom.z() * rotateMatrix[0][2],
                                atom.x() * rotateMatrix[1][0] + atom.y() * rotateMatrix[1][1] + atom.z() * rotateMatrix[1][2],
                                atom.x() * rotateMatrix[2][0] + atom.y() * rotateMatrix[2][1] + atom.z() * rotateMatrix[2][2]
                        );
                    }
                }
                for (Atom atom : residue.getAtoms()) {
                    if (atom.nowhere()) {
                        String name = atom.name;
                        Atom templateAtom = null;
                        for (int i = 0; i < templateResidue.getAtoms().size() & (templateAtom == null); i++) {
                            if (templateResidue.getAtoms().get(i).name().equals(name))
                                templateAtom = templateResidue.getAtoms().get(i);
                        }
                        if ((templateAtom != null) && (!templateAtom.nowhere())) {
                                if (Double.isNaN(templateAtom.x()))
                                    throw new RuntimeException("This is weird.");
                                atom.setXYZ(templateAtom.x(), templateAtom.y(), templateAtom.z());
                            }
                    }
                }
                for (Atom atom : residue.getAtoms()) {
                    if ((atom != null) & (!atom.nowhere()))
                        atom.setXYZ(atom.x() + residueCOM.x(), atom.y() + residueCOM.y(), atom.z() + residueCOM.z());
                }
            }
        }
    }

    private static AtomList getBBAtoms(Residue residue) {
        Atom residueCA  = residue.ca();
        Atom residueCB  = residue.cb();
        Atom residueN   = residue.amideN();
        Atom residueC   = residue.carbonylC();
        AtomList residueAtomsList = new AtomList(residueCA.molecularSystem);
        residueAtomsList.add(residueCA);
        residueAtomsList.add(residueCB);
        residueAtomsList.add(residueN);
        residueAtomsList.add(residueC);
        return residueAtomsList;
    }
    private static Residue getTemplateResidue(Protein template, ResidueType type) {
        for (Residue residue : template.residues())
            if (residue.type() == type)
                return residue;
        throw new RuntimeException("This is weird.");
    }
}

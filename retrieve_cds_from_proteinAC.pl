#!/usr/bin/perl -w
#
# retrieve_cds_from_proteinAC.pl
#
# Description:  a Perl script that takes a list file containing protein accessions from 
#               NCBI (GI or accession number) separated by newlines and retrieves their 
#               CDS sequences. All retrieved sequences are written in "sequence_cds.fna".
#
# Dependencies: Perl;
#               Perl modules listed below (See this link for instructions to install Perl 
#               modules: http://www.cpan.org/modules/INSTALL.html):
#               - Net::SSLeay;
#               - HTTP::Tiny;
#               - XML::Simple.
#
# Author:       Tetsu Sakamoto
#
# Affiliation:  Universidade Federal de Minas Gerais (Brazil)
#               Laboratório de Biodados (ICB, N4-202)
#
# Contact:      tetsufmbio@gmail.com
#

use strict;
use Net::SSLeay;
use HTTP::Tiny;
use XML::Simple qw(XMLin);
use Data::Dumper;

if (scalar @ARGV != 1){
    print "Usage: perl retrieve_cds_from_proteinAC.pl <accession_list_file>\n\n";
	print "<accession_list_file> is a file containing the protein accession from NCBI (GI or accession number) separated by newlines.\n\n" ;
	exit;
}

my $listFile = shift @ARGV;

open(LIST, "< $listFile") or die "Could not open $listFile...\n\n";

my @allList;
my %hash_ids;
my %mapIds;
# check accessions in the list
print "# Checking accessions from the list.\n";
while(my $line = <LIST>){
    chomp $line;
    $line =~ s/[\r\s]//g;
    next if ($line =~ /^$/);
    if ($line =~ /^\d+$/){
        push(@allList, $line);
        $hash_ids{$line}{"type"} = "gi";
    } elsif ($line =~ m/^(NP|AP|XP|YP|WP|ZP)_\d+(\.\d+)?$/ || $line =~ m/^\w{3}\d{5}(\.\d+)?$/){
        push(@allList, $line);
        $hash_ids{$line}{"type"} = "accession";
    } else {
        print "Could not recognize $line as protein identifier from NCBI. Skipped it.\n";
    }
}

# retrieve GI, accession number and taxonomic data of each identifier in the list

my $n = -1;
my $m = -50;
if (scalar @allList > 0){
    print "\n# Retrieving GI, accession number and taxonomic data.\n";
    do {
        $n = $n + 50;
        $m = $m + 50;
        $n = $#allList if ($n > $#allList);
        
        my $url_fetch_id = "https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=text&rettype=seqid&id=".join(",",@allList[$m .. $n]);
        my $fetch_lineage2;
        my $errorCount2 = -1;
        do {
            my $response = HTTP::Tiny->new->get($url_fetch_id);
            $fetch_lineage2 = $response->{content};
            $errorCount2++;
            sleep 1;
        } while ($fetch_lineage2 =~ m/<\/ERROR>|<\/Error>|<title>Bad Gateway!<\/title>|<title>Service unavailable!<\/title>|Error occurred:/ and $errorCount2 < 5);
        if ($errorCount2 > 4){
            die "\nERROR: Sorry, access to NCBI server retrieved error 4 times. Please, try to run this script again later.";
        }
        
        my @ids = split(/\n\n/, $fetch_lineage2);
        foreach my $ids (@ids){
            $ids =~ /accession \"([^\" ]+)\" ,/;
            my $acc = $1;
            $ids =~  /version (\d+) /;
            my $ver = $1;
            $ids =~  /Seq-id ::= gi (\d+)/;
            my $gi = $1;
            $mapIds{"full"}{$gi} = $acc.".".$ver;
            $mapIds{"full"}{$acc.".".$ver} = $gi;
            $mapIds{"part"}{$acc}{"versions"}{$ver} = $gi;
            if (exists $mapIds{"part"}{$acc}{"maxversion"}){
                if ($mapIds{"part"}{$acc}{"maxversion"} < $ver){
                    $mapIds{"part"}{$acc}{"maxversion"} = $ver;
                }
            } else {
                $mapIds{"part"}{$acc}{"maxversion"} = $ver;
            }
        }
        
        # pick xml
        my %map_tax;
        my $url_fetch_seq2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&retmode=xml&id=".join("&id=",@allList[$m .. $n]);
        my $response2 = HTTP::Tiny->new->get($url_fetch_seq2);
        my $link_xml2;
        $errorCount2 = -1;
        do {
            my $response2 = HTTP::Tiny->new->get($url_fetch_seq2);
            $link_xml2 = $response2->{content};
            $errorCount2++;
            sleep 1;
        } while ($link_xml2 =~ m/<\/Error>|<title>Bad Gateway!<\/title>|<title>Service unavailable!<\/title>|Error occurred:/ and $errorCount2 < 5);
        if ($errorCount2 > 4){
            die "\nERROR: Sorry, access to NCBI server retrieved error 4 times. Please, try to run this script again later.";
        }
        
        my $xs1 = XML::Simple->new();
        my $doc_link = $xs1->XMLin($link_xml2, ForceArray => [ 'TSeq' ]);
                
        my @linkSet = @{$doc_link->{"TSeq"}};
        foreach my $link(@linkSet){
            my $accFull = $link->{"TSeq_accver"};
            my $accPart = $accFull;
            $accPart =~ s/\.\d+$//;
            my $tax = $link->{"TSeq_orgname"};
            $map_tax{$accFull} = $tax;
            $map_tax{$accPart} = $tax;
        }
        
        
        foreach my $id (@allList[$m .. $n]){
            if (exists $mapIds{"full"}{$id}){
                $hash_ids{$id}{"convert"} = $mapIds{"full"}{$id};
            } elsif (exists $mapIds{"part"}{$id}) {
                my $maxVer = $mapIds{"part"}{$id}{"maxversion"};
                $hash_ids{$id}{"convert"} = $mapIds{"part"}{$id}{"versions"}{$maxVer};
            } else {
                print "Could not retrieve some data from $id. Skipped it.\n";
                delete $hash_ids{$id};
            }
            
            if ($hash_ids{$id}{"type"} eq "accession"){
                if (exists $map_tax{$id}){
                    $hash_ids{$id}{"orgname"} = $map_tax{$id};
                } else {
                    print "Could not retrieve some data from $id. Skipped it.\n";
                    delete $hash_ids{$id};
                }
            } else {
                if (exists $map_tax{$hash_ids{$id}{"convert"}}){
                    $hash_ids{$id}{"orgname"} = $map_tax{$hash_ids{$id}{"convert"}};
                } else {
                    print "Could not retrieve some data from $id. Skipped it.\n";
                    delete $hash_ids{$id};
                }
            }
        }
    } while ($n < $#allList);
}

my @mrna;
my @nomrna;
my @genome;
$n = -1;
$m = -50;

if (scalar @allList > 0){
    print "\n# Linking each protein accession to nuccore database.\n";
    do {
        $n = $n + 50;
        $m = $m + 50;
        $n = $#allList if ($n > $#allList);
        my $url_fetch_seq = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=nuccore&id=".join("&id=",@allList[$m .. $n]);
        my $response = HTTP::Tiny->new->get($url_fetch_seq);
        my $link_xml;
        my $errorCount = -1;
        do {
            my $response = HTTP::Tiny->new->get($url_fetch_seq);
            $link_xml = $response->{content};
            $errorCount++;
            sleep 1;
        } while ($link_xml =~ m/<\/Error>|<title>Bad Gateway!<\/title>|<title>Service unavailable!<\/title>|Error occurred:/ and $errorCount < 5);
        if ($errorCount > 4){
            die "\nERROR: Sorry, access to NCBI server retrieved error 4 times. Please, try to run this script again later.";
        }
        my $xs1 = XML::Simple->new();
        my $doc_link = $xs1->XMLin($link_xml, ForceArray => [ 'LinkSet' ]);

        my %hashTemp;
        my @linkSet = @{$doc_link->{"LinkSet"}};
        foreach my $link(@linkSet){
            my $protID = $link->{"IdList"}->{"Id"};
            my $mrna;
            my $nomrna;
            
            if (exists $link->{"LinkSetDb"}){
                if (ref $link->{"LinkSetDb"} eq 'ARRAY'){
                    my @linkSetDb = @{$link->{"LinkSetDb"}};
                    foreach my $linkSetDb (@linkSetDb){
                        
                        if ($linkSetDb->{"LinkName"} eq "protein_nuccore_mrna"){
                            $mrna = $linkSetDb->{"Link"}->{"Id"};
                            $nomrna = "";
                            last;
                        } elsif ($linkSetDb->{"LinkName"} eq "protein_nuccore"){
                            if (ref $linkSetDb->{"Link"} eq 'ARRAY'){
                                $nomrna = $linkSetDb->{"Link"}->[0]->{"Id"};
                            } else {
                                $nomrna = $linkSetDb->{"Link"}->{"Id"};
                            }
                        }
                    }
                } else {
                    if ($link->{"LinkSetDb"}->{"LinkName"} eq "protein_nuccore_mrna"){
                        $mrna = $link->{"LinkSetDb"}->{"Link"}->{"Id"};                            
                    } elsif ($link->{"LinkSetDb"}->{"LinkName"} eq "protein_nuccore"){
                        if (ref $link->{"LinkSetDb"}->{"Link"} eq 'ARRAY'){
                            $nomrna = $link->{"LinkSetDb"}->{"Link"}->[0]->{"Id"};
                        } else {
                            $nomrna = $link->{"LinkSetDb"}->{"Link"}->{"Id"};
                        }
                    }
                }
            }
            if ($mrna){
                push(@mrna, $mrna);
                $hashTemp{$protID}{"mrna"} = $mrna;
            } elsif ($nomrna) {
                $hashTemp{$protID}{"genome"} = $nomrna;
                push(@genome, $nomrna);
                push(@nomrna, $protID);
            }
        }
        
        foreach my $id (@allList[$m .. $n]){
            my $searchID = $id;
            if ($hash_ids{$id}{"type"} eq "accession"){
                $searchID = $hash_ids{$id}{"convert"};
            }
            
            if (exists $hashTemp{$searchID}){
                if (exists $hashTemp{$searchID}{"mrna"}){
                    $hash_ids{$id}{"mrna"} = $hashTemp{$searchID}{"mrna"};
                } else {
                    $hash_ids{$id}{"genome"} = $hashTemp{$searchID}{"genome"};
                }
            } else {
                print "Could not retrieve nucleotide data from $id. Skipped it.\n";
                delete $hash_ids{$id};
            }
        }
    }  while ($n < $#allList);
}

@allList = keys %hash_ids;
my %hash_mrna;
$n = -1;
$m = -50;

if (scalar @mrna > 0){
    print "\n# Retrieving mrna sequences.\n";
    do {
        $n = $n + 50;
        $m = $m + 50;
        $n = $#mrna if ($n > $#mrna);

        # pick ids
        my $url_fetch_id = "https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&retmode=text&rettype=seqid&id=".join(",",@mrna[$m .. $n]);
        my $fetch_lineage2;
        my $errorCount2 = -1;
        do {
            my $response = HTTP::Tiny->new->get($url_fetch_id);
            $fetch_lineage2 = $response->{content};
            $errorCount2++;
            sleep 1;
        } while ($fetch_lineage2 =~ m/<\/ERROR>|<\/Error>|<title>Bad Gateway!<\/title>|<title>Service unavailable!<\/title>|Error occurred:/ and $errorCount2 < 5);
        if ($errorCount2 > 4){
            die "\nERROR: Sorry, access to NCBI server retrieved error 4 times. Please, try to run this script again later.";
        }
        
        my @ids = split(/\n\n/, $fetch_lineage2);
        foreach my $ids (@ids){
            $ids =~ /accession \"([^\" ]+)\" ,/;
            my $acc = $1;
            $ids =~  /version (\d+) /;
            my $ver = $1;
            $ids =~  /Seq-id ::= gi (\d+)/;
            my $gi = $1;
            $mapIds{"full"}{$gi} = $acc.".".$ver;
            $mapIds{"full"}{$acc.".".$ver} = $gi;
            $mapIds{"part"}{$acc}{"versions"}{$ver} = $gi;
            if (exists $mapIds{"part"}{$acc}{"maxversion"}){
                if ($mapIds{"part"}{$acc}{"maxversion"} < $ver){
                    $mapIds{"part"}{$acc}{"maxversion"} = $ver;
                }
            } else {
                $mapIds{"part"}{$acc}{"maxversion"} = $ver;
            }
        }
        
        # pick fasta
        my $url_fetch_seq = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta_cds_na&retmode=text&id=".join("&id=",@mrna[$m .. $n]);
        my $response = HTTP::Tiny->new->get($url_fetch_seq);
        my $link_xml;
        my $errorCount = -1;
        do {
            my $response = HTTP::Tiny->new->get($url_fetch_seq);
            $link_xml = $response->{content};
            $errorCount++;
            sleep 1;
        } while ($link_xml =~ m/<\/Error>|<title>Bad Gateway!<\/title>|<title>Service unavailable!<\/title>|Error occurred:/ and $errorCount < 5);
        if ($errorCount > 4){
            die "\nERROR: Sorry, access to NCBI server retrieved error 4 times. Please, try to run this script again later.";
        }
        my @fasta = split(/^>|[\n]+>/, $link_xml);
        my %hashTemp;
        for(my $i = 1; $i < scalar @fasta; $i++){
            my ($header, $seq) = split(/\n/, $fasta[$i], 2);
            $header =~ /lcl\|([\w\.]+)_cds_/;
            my $acc = $1;
            $hashTemp{$acc} = $seq;
            #print OUT ">".$map_tax{$acc}."\n".$seq;
        }
        foreach my $mrna(@mrna[$m .. $n]){
            if (exists $mapIds{"full"}{$mrna}){
                my $defAcc = $mapIds{"full"}{$mrna};
                if (exists $hashTemp{$defAcc}){
                    $hash_mrna{$mrna} = $hashTemp{$defAcc};
                }
            }
        }
    } while ($n < $#mrna);
    
}

$n = -1;
$m = -50;
my %hash_nomrna;
if (scalar @nomrna > 0){
    print "\n# Retrieving cds from genome data.\n";
    my %accessions;
    
    foreach my $nomrna(@nomrna){
        $accessions{$mapIds{"full"}{$nomrna}} = $nomrna;
    }
    
    $n = -1;
    $m = -25;
    do {
        $n = $n + 25;
        $m = $m + 25;
        $n = $#genome if ($n > $#genome);
        
        # pick accession
        my %map_tax;
        my $url_fetch_seq2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta_cds_na&retmode=text&id=".join("&id=",@genome[$m .. $n]);
        my $response2 = HTTP::Tiny->new->get($url_fetch_seq2);
        my $link_xml2;
        my $errorCount2 = -1;
        do {
            my $response2 = HTTP::Tiny->new->get($url_fetch_seq2);
            $link_xml2 = $response2->{content};
            $errorCount2++;
            sleep 1;
        } while ($link_xml2 =~ m/<\/Error>|<title>Bad Gateway!<\/title>|<title>Service unavailable!<\/title>|Error occurred:/ and $errorCount2 < 5);
        if ($errorCount2 > 4){
            die "\nERROR: Sorry, access to NCBI server retrieved error 4 times. Please, try to run this script again later.";
        }
        
        my @link = split(">", $link_xml2);
        my $discard = shift @link;
        foreach my $fasta(@link){
            my ($header, $seq) = split(/\n/, $fasta, 2);
            $header =~ /\[protein_id=(\wP_\d+\.\d+)\]/;
            my $id = $1;
            if (exists $accessions{$id}){
                #print OUT ">".$header."\n".$seq;
                $hash_nomrna{$accessions{$id}} = $seq;
            }
        }
        
    } while ($n < $#genome);
}

if (scalar @allList > 0){
    open(OUT, "> sequence_cds.fna") or die "Could not create file sequence_cds.fna.\n\n";
    @allList = keys %hash_ids;
    foreach my $id (@allList){
        if (exists $hash_ids{$id}{"mrna"}){
            my $mrna = $hash_ids{$id}{"mrna"};
            if (exists $hash_mrna{$mrna}){
                print OUT ">".$id." ".$hash_ids{$id}{"orgname"}."\n".$hash_mrna{$mrna};
            } else {
                print "Could not retrieve mrna sequence from $id.\n";
            }
        } elsif (exists $hash_ids{$id}{"genome"}){
            my $searchID = $id;
            if ($hash_ids{$id}{"type"} eq "accession"){
                $searchID = $hash_ids{$id}{"convert"};
            }
            if (exists $hash_nomrna{$searchID}){
                print OUT ">".$id." ".$hash_ids{$id}{"orgname"}."\n".$hash_nomrna{$searchID};
            } else {
                print "Could not retrieve cds sequence from $id.\n";
            }
        } else {
            print "Could not retrieve nucleotide sequence from $id.\n";
        }
    }
    close OUT;
    print "\n# CDS sequence printed in sequence_cds.fna file.\n";
}

print "\n# All done!\n";
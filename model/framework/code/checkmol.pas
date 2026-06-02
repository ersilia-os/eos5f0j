program checkmol;
(*
checkmol/matchmol
Norbert Haider, University of Vienna, 2003-2018
norbert.haider@univie.ac.at

This software is published under the terms of the GNU General Public
License (GPL, see below). For a detailed description of this license,
see http://www.gnu.org/copyleft/gpl.html

If invoked as "checkmol", this program reads 2D and 3D molecular structure
files in various formats (Tripos Alchemy "mol", Tripos SYBYL "mol2", MDL "mol")
and describes the molecule either by its functional groups or by a set of
descriptors useful for database pre-screening (number of rings, sp2-hybridized
carbons, aromatic bonds, etc.; see definition of molstat_rec).

If invoked as "matchmol", the program reads two individual 2D or 3D molecular
structure files in various formats (see above) and checks if the first molecule
(the "needle") is a substructure of the second one (the "haystack").
"Haystack" can also be a MDL SD-file (containing multiple MOL files);
if invoked with "-" as file argument, both "needle" and "haystack" are
read as an SD-file from standard input, assuming the first entry in
the SDF to be the "needle"; output: entry number + ":F" (false) or ":T" (true)


Compile with fpc (Free Pascal, see http://www.freepascal.org), using
the -Sd or -S2 option (Delphi mode; IMPORTANT!)

example for compilation (with optimization) and installation:

fpc checkmol.pas -S2 -O3 -Op3

as "root", do the following:

cp checkmol /usr/local/bin
cd /usr/local/bin
ln checkmol matchmol         (ATTENTION: a symbolic link does not work!)


Version history

v0.1   extends "chkmol" utility (N. Haider, 2002), adds matching functionality;

v0.1a  minor bugfixes:
       stop ring search when max_rings is reached (function is_ringpath);
       fixed upper/lowercase of element symbol in MDL molfile output;
       added debug output for checkmol (-D option)

v0.2   added functionality to switch ring search from SAR (set of all rings) to
       SSR (set of small rings; see below) if number of rings exceeds max_rings
       (1024), e.g. with Buckminster-fullerenes (thanks to E.-G. Schmid for a hint);
       added -r command-line option to force ring search into SSR mode;
       as SSR, we regard the set of all rings with max. 10 ring members, and in
       which no ring is completely contained in another one

v0.2a  fixed a bug in function is_ringpath which could cause SAR ring search to
       overlook a few rings (e.g., the six 8-membered rings in cubane)

v0.2b  fixed unequal treatment of needle and haystack structures with respect
       to aromaticity check (trusting input file vs. full check)

v0.2c  modified the changes of v0.2b so that they affect only "matchmol" mode;
       added quick_match function

v0.2d  refined function bondtypes_OK so that it does not accept C=O fragments
       (and C=S, C=Se, C=Te) as equivalent to aromatic C-O (C-S, C-Se, C-Te);
       fixed function find_ndl_ref_atom to use only heavy atoms as "needle"
       reference atoms for atom-to-atom match; update function quick_match to
       handle queries with only one heavy atom (thanks to A. Barozza for a hint)

v0.2e  modified procedure get_molstat: adds 1 formal double bond for each aromatic 
       ring to the "fingerprint" in order to allow an isolated double bond in the 
       needle to be matched as a fragment of an aromatic ring
       ATTENTION: "fingerprints" (molecular statistics) generated with a previous
       version of checkmol(-x and -X option) must be updated for structure/
       substructure database searching with matchmaol.

v0.2f  modified procedures chk_ccx, chk_xccx, write_fg_text, write_fg_text_de,
       write_fg_code in order to report halogen derivatives other than alkyl halides,
       aryl halides and acyl halides (e.g., vinyl halides) and C=C compounds other
       than "true" alkenes, enols, enediols, enamines, etc. (e.g. vinyl halides);
       added "strict mode" (option -s) for matching: this uses a more thorough
       comparison of atom types and bond types (new functions atomtypes_OK_strict,
       bondtypes_OK_strict, several minor modifications in other subroutines).

v0.2g  changed procedure readinputfile (+ minor changes in main routine) in
       order to read very large SD input files without "not enough memory" error
       (the previous version attempted to read the entire SD file into a buffer,
       the new version reads only one molecule at once);
       fixed a minor bug in read_mol2file which led to undefined bond types;
       added definition of "debug" option as a compiler flag.

v0.2h  fixed a small bug which caused program hangs with (e.g.) some metal 
       complexes: this was solved just by increasing the number of possible 
       neighbor atoms from 8 to 16 (now defined in the constant max_neighbors).       

v0.2i  introduced some more plausibility checking (number of directly connected 
       atoms), result stored in variable mol_OK (set in count_neighbors); 
       completely rewrote function matrix_OK which now can handle up to 
       max_neighbors substituents per atom (previous versions were limited to
       4 substituents; larger numbers are required e.g. for ferrocenes and 
       similar compounds).

v0.2j  new function raw_hetbond_count: ignores bond order (in contrast to function
       hetbond_count), this is better suitable for molecular statistics.
       ATTENTION: previously generated "fingerprints" in databases must be updated!
       improved recognition of non-conformant SD files (with missing "M  END" line)

v0.2k  changed quick_match routine to compare atom types only by element in
       order to avoid (rare) cases of non-matching identical input files;
       slightly changed error message for non-existant input files

v0.3   changed function update_Htotal in order to distinguish between 3-valent
       and 5-valent phosphorus (thanks to H. Feldman for this suggestion);
       added a table (array ringprop) to store ring sizes and aromaticity for
       faster lookup; changed aromaticity detection (chk_arom) to be fully
       independent of Kekulé structures in condensed ring systems; changed add_ring
       to store new rings in ascending order (with respect to ring size): this
       will cause the aromaticity detection to start with smaller rings;
       added additional calls to chk_arom when in SSR ring search mode (to ensure
       that all aromatic rings are found); speeded up function is_newring;
       added option "-M": this restores the behavior of previous versions,
       i.e. metal atoms will be accepted as ring members (which costs a lot
       of performance for coordinate compounds like ferrocenes but gives
       only little useful information); when used in databases: use the _same_
       mode ("-M" or not) _both_ for checkmol (creation of fingerprints) and matchmol;
       fixed ugly linebreaks in show_usage;

v0.3a  fixed a bug in read_charges (which was introduced in the v0.3 code cleanup)

v0.3b  fixed recognition of ketenes, CO2, CS2, CSO, phosphines, boronic esters

v0.3c  fixed recognition of hydrazines, hydroxylamines, thiocarboxylic acids,
       orthocarboxylic acids; slightly changed textual output for a few functional
       groups

v0.3d  added bond topology feature ("any", "ring", "chain", as defined by MDL specs, 
       plus "always_any", "excess-ringcount", "exact-ringcount") to bond properties 
       and implemented its use for substructure searches, either as specified in the 
       query MDL molfile or by using "strict" mode; added ring_count property to 
       "bonds" section of Checkmol-tweaked molfiles (using an unused field in the MDL 
       molfile specs);
       added option for E/Z isomer search, either globally (-g option) or per
       double bond: bstereo_xyz property, encoded either by a leading "0" in the 
       bond block of a checkmol-tweaked molfile or by using the "stereo care" flag
       as defined in the MDL file specs (both atoms of a double bond must carry
       this flag)

v0.3e  fixed wrong position of "stereo care" flag expected in input molfile; 

v0.3f  added option for enantiospecific search (R/S isomer search); this can be
       turned on either globally (-G option or using the "chiral" flag in the
       "counts line" of the query MDL molfile) or per atom (using - or abusing -
       the "stereo care" flag); enantiomers are compared using the 3D coordinates
       of the substituents at a chiral center, so we do not have to rely on the
       presence of correct "up" and "down" bond types in the candidate structures;
       nevertheless, "pseudo-3D" structures (i.e. flat 2D structures with "up"
       and "down" bond notation) can be also used, even in mixed mode both as
       query structure and candidate structure; added support for "up" and "down"
       bond notation in MDL molfile output (if -m option is used and these bond 
       types were present in the input file); refined function find_ndl_ref_atom;
       fixed Alchemy atom type misinterpretation (upper/lower case mismatch);

v0.3g  minor fixes: recognition of sp2 nitrogen in N=N structures (function 
       update_atypes); extended E/Z geometry check to C=N and N=N bonds; accelerated
       exact search by initial comparison of C,O,N counts; added meaningful error
       message for input structures with 0 (zero) atoms;

v0.3h  added a match matrix clean-up step to function is_matching in order to
       remove "impossible" multiple matches prior to chirality check; added
       function ndl_maybe_chiral as a plausibility check; added support for other
       atoms than carbon as chiral centers; changed function is_cis from a
       distance-based cis/trans check into a torsion-based check;

v0.3i  minor fixes in functions max_neighbors and chk_ammon; 
       revert ringsearch_mode to its original value for each new molecule 
       (main routine)

v0.3j  various improvements in checkmol: added alkenyl (vinyl) and alkynyl
       residues as possible substituents for amides, esters, ethers, etc.; 
       extended recognition of orthocarboxylic acid derivatives; refined 
       recognition of oxohetarenes and iminohetareenes, ureas, nitroso 
       compounds; fixed reading of "old-style" charges from MDL files; 
       refined aromaticity check; fixed a bug in procedure chk_arom and 
       renamed function is_exocyclic_methylene_C into 
       find_exocyclic_methylene_C; added assignment of pointers to "nil" 
       after calling "freemem" in procedure zap_molecule and zap_needle;
       matchmol: refined selection of needle reference atom; added "strict 
       chirality check" mode (if -x -s and -G options are used): returns 
       "false" if a chiral or pseudochiral atom is matched against its 
       counterpart with undefined or "flat" geometry; added an alternative 
       selection mechanism for the needle reference atom, based on the 
       Morgan algorithm;

v0.3k  improved matching of nitro compounds (ionic vs. non-ionic 
       representation); some finetuning in recognition of N-oxides,
       hydroxylamines, etc.; modified get_molstat in order to ignore 
       charges in nitro groups, N-oxides, sulfoxides, azides, etc. in 
       ionic representation; added function normalize_ionic_bonds; 
       fixed a bug in is_alkenyl; fixed a bug (introduced in v0.3j) in 
       the treatment of "old-style" MDL charges;
       
v0.3l  minor adjustments in quick_match (for unknown atom types,
       opposite charges); some performance improvements, e.g. in function 
       path_pos, reduced extensive calls to is_heavyatom and is_metal 
       by extending atom_rec with boolean fields 'heavy' and 'metal' 
       (thanks to E.-G.Schmid); added molstat field n_rBz (number of 
       benzene rings; deactivated by default); added -l command-line 
       option to checkmol (lists all molstat codes with a short 
       explanation); added graceful handling of NoStruct (i.e.,
       zero-atom) molecules

v0.3m  minor bug fixes (chk_imine), let normalize_ionic_bonds return "true"
       if some bond-order change was made; fixed incorrect implementation 
       of the Morgan algorithm (thanks to H.Feldman), fixed is_matching in 
       order to prevent wrong matches of larger rings on smaller ones; count 
       halogens also as type "X" in molstat (as before v0.3l); added some more
       molstat descriptors for most elements (deactivated by default; for 
       activation, uncomment the definition of "extended_molstat", see below);
       added support for the generation of binary fingerprints with output
       either as boolean values (-f option) or as decimal values (-F option), 
       fingerprint patterns are supplied by the user (in SDF format, max. 62 
       structures per file); added a version check ("tweaklevel", 
       ringsearch_mode, opt_metalrings) for "tweaked" MDL molfiles
       
v0.3n  increased max_vringsize in SSR mode from 10 to 12 (now defined in
       ssr_vringsize); optionally changed ring statistics to ignore envelope 
       rings (which completely contain another ring, "reduced SAR") in order to 
       make molstat in SAR and SSR mode more compatible) - disabled by default; 
       introduced a new molstat descriptor: n_br2p (number of bonds which 
       belong to two or more rings); changed n_b1 counter to ignore bonds to 
       explicit hydrogens; added procedure fix_ssr_ringcounts
       (see comment in the code); added bond property mdl_stereo for preservation
       of original value; slightly changed "get_molstat" and "update_atypes" in 
       order to consolidate atom type for various N atoms; remember "chiral" flag
       status in molfile output (-m); modified chk_arom in order to accept rings
       as aromatic if all bonds are of type 'A'; several minor bug fixes

v0.3o  minor changes in update_atypes (nitrogen); changed the matching 
       algorithm in order to take care of disconnected molecular graphs 
       (e.g., salts); refined matching of atoms with formal charges
       
v0.3p  added a maximum recursion depth in function is_ringpath() as another 
       SAR->SSR fallback criterion in order to speed up execution time for 
       molecules with _many_ fused rings (thanks to E.-G. Schmid for this 
       suggestion); added the same recursion depth lilitation to is_matching();
       refined charge comparison and added support for isotope and radical type 
       matching (this is actually a back-port from E.-G. Schmid's "Barsoi", 
       the C version); some other minor fixes, including better matching for
       rings containing wildcard atoms (Q, A) or bonds (MDL type 5, 6, 7, 8)

v0.4   added support for hash-based fingerprints (defaults: 512 bits, 2 bits
       per fragment, minimum fragment length = 3 atoms, maximum fragment length =
       8 atoms); changed output format for -b option into the decimal representation
       (32-bit unsigned integers) of the bitstring representing the functional
       groups; allowed more output options to be combined; lots of minor bugfixes       

v0.4a  added a global match matrix + extended output option (-n, -N) for matchmol: 
       hits are indicated as (e.g.) 1:T 1=3;2=5,6;3=6,5 which means: structure no.1
       is a hit (:T) and needle atom 1 matches haystack atom 3, needle atom 2
       matches haystack atoms 5 and 6, needle atom 3 again matches haystack 
       atoms 6 and 5; -n shows all matching atom pairs, -N only those for
       the first matching substructure which is found in the candidate structure;
       minor bug fix in checkmol: corrected recognition of primary thiocarboxamides;

v0.4b  added procedure update_atypes_quick; fixed bondtypes_OK(), bondtypes_OK_strict()
       and assemble_frag() in order to distinguish between arenes and dehydroarenes 
       (i.e., aromatic rings with a triple bond); added a count for true heavy atoms
       (ignoring deuterium and tritium) to avoid false positive hits in -x mode; other
       minor bug fixes, including a fix to avoid overlapping matches (in -s, -n and -N 
       mode only) and an improvement of is_alkenyl(); changed output of -n option in 
       order to avoid loss of information; added support for a separator label 
       (comment) in a SDF $$$$ line

v0.4c  fixed a minor bug in find_ndl_ref_atom which was introduced in v0.4b;
       fixed a minor bug in chk_c_n(); minor changes in handling of sep_label;

v0.4d  added "xx" option for full molstat output with field names; added calculation
       of molecular formula and molecular weight (output in line 3 of MDL molfile);
       minor bug fix in find_ndl_ref_atom;
       
v0.4e  fixed a bug in command-line parsing; fixed a collision in the hashed
       fingerprints of halogen-containing fragments; added "nochirality" option; 
       added support for some branched fragments (C and N atoms only) in the
       hashed fingerprints; ATTENTION: please re-generate any hashed-fingerprint 
       tables in MolDB5R if you upgrade to this version of checkmol/matchmol!

v0.5   changed the internal handling of clear-text functional group names in order 
       to be more flexible towards additional languages; added a new option (-p)
       which lists the number and position of all detected functional groups in
       a molecule (for most groups, a single key atom is listed, for some groups
       a connected pair of key atoms is used, for aromatic and heterocyclic rings
       all ring atoms are listed; for a detailed description, see the document
       fgtable.pdf, available from the checkmol/matchmol homepage)

v0.5a  fixed an issue with charge handling in read_mol2file(); refined function
       normalize_ionic_bonds to increase bond order from double to triple to
       only C=C and C=N bonds; minor change in atomtypes_OK_strict()

v0.5b  fixed an issue with query bond type "any" in JSME (which uses "9" instead 
       of "8"); minor changes in chk_c_n(), now recognizing nitro groups etc. attached
       to any kind of carbon

Credits:

Rami Jbara (rami_jbara@hotmail.com) designed the 8-character functional
group codes.

I am also very grateful to Ernst-Georg Schmid (Bayer Business Services, 
Leverkusen/Germany), to Howard Feldman (The Blueprint Initiative, 
Toronto/Canada), and to Robert Bruccoleri (Congenomics) for numerous bug 
reports and suggestions for improvements.



===============================================================================
DISCLAIMER
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
===============================================================================
*)

//{$DEFINE debug}                    // uncomment if desired
//{$DEFINE extended_molstat}         // uncomment if desired
//{$DEFINE reduced_SAR}              // uncomment if desired
{$DEFINE extended_hfp}             // uncomment if desired (hfp for some branched fragments)

{$NOTES OFF}
{$WARNINGS OFF}

uses
  SYSUTILS, MATH;


const
  version       = '0.5b';
  tweaklevel    = 2;   // v0.3m; accept tweaks in MDL molfiles only if same level
  max_atoms     = 1024;
  max_bonds     = 1024;
  max_ringsize  = 128;
  max_rings     = 1024;
  max_fg        = 256; // must be a multiple of 32
  used_fg       = 204; // v0.5; number of functional groups currently used (204 in v0.5)
  max_neighbors = 20;  // new in v0.2h (16); increased in v0.4d (20)
  TAB           = #26;

  max_matchpath_length = 256;
  pmCheckMol   = 1001;
  pmMatchMol   = 1002;

  rs_sar       = 2001;   // ring search mode: SAR = set of all rings
  rs_ssr       = 2002;   //                   SSR = set of small rings

  btopo_any        = 0;  // bond topology, new in v0.3d
  btopo_ring       = 1;  //
  btopo_chain      = 2;  //
  btopo_always_any = 3;  // even in "strict mode"
  btopo_excess_rc  = 4;  // bond in candidate must be included in _more_ rings than
                         // the matching bond in the query ==> specific search for
                         // annulated systems
  btopo_exact_rc   = 5;  // bond in query and candidate must have same ring count

  bstereo_any      = 0;  // new in v0.3d, any E/Z isomer matches (for double bonds)
  bstereo_xyz      = 1;  // E/Z match is checked by using XYZ coordinates of the atoms
  bstereo_up       = 11; // new in v0.3f, flags for single bonds
  bstereo_down     = 16; //
  fpf_boolean      = 3001; // v0.3m; format for fingerprint output (n:T, n:F)
  fpf_decimal      = 3002; // v0.3m; format for fingerprint output as decimal value of bitstring
  fp_blocksize     = 62;   // v0.3m; 1 bit may be used for sign (e.g. by PHP), 1 bit for "exact hit"
  ssr_vringsize    = 12;   // v0.3n; max. ring size in SSR mode
  max_recursion_depth       = 500000;  // v0.3p
  max_match_recursion_depth = 500000;  // v0.3p
  max_ndl_gmmsize  = 128;  // v0.4a
  max_hst_gmmsize  = 1024; // v0.4a


  // Definitions for functional groups:
  fg_cation                         = 001;
  fg_anion                          = 002;
  fg_carbonyl                       = 003;
  fg_aldehyde                       = 004;
  fg_ketone                         = 005;
  fg_thiocarbonyl                   = 006;
  fg_thioaldehyde                   = 007;
  fg_thioketone                     = 008;
  fg_imine                          = 009;
  fg_hydrazone                      = 010;
  fg_semicarbazone                  = 011;
  fg_thiosemicarbazone              = 012;
  fg_oxime                          = 013;
  fg_oxime_ether                    = 014;
  fg_ketene                         = 015;
  fg_ketene_acetal_deriv            = 016;
  fg_carbonyl_hydrate               = 017;
  fg_hemiacetal                     = 018;
  fg_acetal                         = 019;
  fg_hemiaminal                     = 020;
  fg_aminal                         = 021;
  fg_thiohemiaminal                 = 022;
  fg_thioacetal                     = 023;
  fg_enamine                        = 024;
  fg_enol                           = 025;
  fg_enolether                      = 026;
  fg_hydroxy                        = 027;
  fg_alcohol                        = 028;
  fg_prim_alcohol                   = 029;
  fg_sec_alcohol                    = 030;
  fg_tert_alcohol                   = 031;
  fg_1_2_diol                       = 032;
  fg_1_2_aminoalcohol               = 033;
  fg_phenol                         = 034;
  fg_1_2_diphenol                   = 035;
  fg_enediol                        = 036;
  fg_ether                          = 037;
  fg_dialkylether                   = 038;
  fg_alkylarylether                 = 039;
  fg_diarylether                    = 040;
  fg_thioether                      = 041;
  fg_disulfide                      = 042;
  fg_peroxide                       = 043;
  fg_hydroperoxide                  = 044;
  fg_hydrazine                      = 045;
  fg_hydroxylamine                  = 046;
  fg_amine                          = 047;
  fg_prim_amine                     = 048;
  fg_prim_aliph_amine               = 049;
  fg_prim_arom_amine                = 050;
  fg_sec_amine                      = 051;
  fg_sec_aliph_amine                = 052;
  fg_sec_mixed_amine                = 053;
  fg_sec_arom_amine                 = 054;
  fg_tert_amine                     = 055;
  fg_tert_aliph_amine               = 056;
  fg_tert_mixed_amine               = 057;
  fg_tert_arom_amine                = 058;
  fg_quart_ammonium                 = 059;
  fg_n_oxide                        = 060;
  fg_halogen_deriv                  = 061;
  fg_alkyl_halide                   = 062;
  fg_alkyl_fluoride                 = 063;
  fg_alkyl_chloride                 = 064;
  fg_alkyl_bromide                  = 065;
  fg_alkyl_iodide                   = 066;
  fg_aryl_halide                    = 067;
  fg_aryl_fluoride                  = 068;
  fg_aryl_chloride                  = 069;
  fg_aryl_bromide                   = 070;
  fg_aryl_iodide                    = 071;
  fg_organometallic                 = 072;
  fg_organolithium                  = 073;
  fg_organomagnesium                = 074;
  fg_carboxylic_acid_deriv          = 075;
  fg_carboxylic_acid                = 076;
  fg_carboxylic_acid_salt           = 077;
  fg_carboxylic_acid_ester          = 078;
  fg_lactone                        = 079;
  fg_carboxylic_acid_amide          = 080;
  fg_carboxylic_acid_prim_amide     = 081;
  fg_carboxylic_acid_sec_amide      = 082;
  fg_carboxylic_acid_tert_amide     = 083;
  fg_lactam                         = 084;
  fg_carboxylic_acid_hydrazide      = 085;
  fg_carboxylic_acid_azide          = 086;
  fg_hydroxamic_acid                = 087;
  fg_carboxylic_acid_amidine        = 088;
  fg_carboxylic_acid_amidrazone     = 089;
  fg_nitrile                        = 090;
  fg_acyl_halide                    = 091;
  fg_acyl_fluoride                  = 092;
  fg_acyl_chloride                  = 093;
  fg_acyl_bromide                   = 094;
  fg_acyl_iodide                    = 095;
  fg_acyl_cyanide                   = 096;
  fg_imido_ester                    = 097;
  fg_imidoyl_halide                 = 098;
  fg_thiocarboxylic_acid_deriv      = 099;
  fg_thiocarboxylic_acid            = 100;
  fg_thiocarboxylic_acid_ester      = 101;
  fg_thiolactone                    = 102;
  fg_thiocarboxylic_acid_amide      = 103;
  fg_thiolactam                     = 104;
  fg_imido_thioester                = 105;
  fg_oxohetarene                    = 106;
  fg_thioxohetarene                 = 107;
  fg_iminohetarene                  = 108;
  fg_orthocarboxylic_acid_deriv     = 109;
  fg_carboxylic_acid_orthoester     = 110;
  fg_carboxylic_acid_amide_acetal   = 111;
  fg_carboxylic_acid_anhydride      = 112;
  fg_carboxylic_acid_imide          = 113;
  fg_carboxylic_acid_unsubst_imide  = 114;
  fg_carboxylic_acid_subst_imide    = 115;
  fg_co2_deriv                      = 116;
  fg_carbonic_acid_deriv            = 117;
  fg_carbonic_acid_monoester        = 118;
  fg_carbonic_acid_diester          = 119;
  fg_carbonic_acid_ester_halide     = 120;
  fg_thiocarbonic_acid_deriv        = 121;
  fg_thiocarbonic_acid_monoester    = 122;
  fg_thiocarbonic_acid_diester      = 123;
  fg_thiocarbonic_acid_ester_halide = 124;
  fg_carbamic_acid_deriv            = 125;
  fg_carbamic_acid                  = 126;
  fg_carbamic_acid_ester            = 127;
  fg_carbamic_acid_halide           = 128;
  fg_thiocarbamic_acid_deriv        = 129;
  fg_thiocarbamic_acid              = 130;
  fg_thiocarbamic_acid_ester        = 131;
  fg_thiocarbamic_acid_halide       = 132;
  fg_urea                           = 133;
  fg_isourea                        = 134;
  fg_thiourea                       = 135;
  fg_isothiourea                    = 136;
  fg_guanidine                      = 137;
  fg_semicarbazide                  = 138;
  fg_thiosemicarbazide              = 139;
  fg_azide                          = 140;
  fg_azo_compound                   = 141;
  fg_diazonium_salt                 = 142;
  fg_isonitrile                     = 143;
  fg_cyanate                        = 144;
  fg_isocyanate                     = 145;
  fg_thiocyanate                    = 146;
  fg_isothiocyanate                 = 147;
  fg_carbodiimide                   = 148;
  fg_nitroso_compound               = 149;
  fg_nitro_compound                 = 150;
  fg_nitrite                        = 151;
  fg_nitrate                        = 152;
  fg_sulfuric_acid_deriv            = 153;
  fg_sulfuric_acid                  = 154;
  fg_sulfuric_acid_monoester        = 155;
  fg_sulfuric_acid_diester          = 156;
  fg_sulfuric_acid_amide_ester      = 157;
  fg_sulfuric_acid_amide            = 158;
  fg_sulfuric_acid_diamide          = 159;
  fg_sulfuryl_halide                = 160;
  fg_sulfonic_acid_deriv            = 161;
  fg_sulfonic_acid                  = 162;
  fg_sulfonic_acid_ester            = 163;
  fg_sulfonamide                    = 164;
  fg_sulfonyl_halide                = 165;
  fg_sulfone                        = 166;
  fg_sulfoxide                      = 167;
  fg_sulfinic_acid_deriv            = 168;
  fg_sulfinic_acid                  = 169;
  fg_sulfinic_acid_ester            = 170;
  fg_sulfinic_acid_halide           = 171;
  fg_sulfinic_acid_amide            = 172;
  fg_sulfenic_acid_deriv            = 173;
  fg_sulfenic_acid                  = 174;
  fg_sulfenic_acid_ester            = 175;
  fg_sulfenic_acid_halide           = 176;
  fg_sulfenic_acid_amide            = 177;
  fg_thiol                          = 178;
  fg_alkylthiol                     = 179;
  fg_arylthiol                      = 180;
  fg_phosphoric_acid_deriv          = 181;
  fg_phosphoric_acid                = 182;
  fg_phosphoric_acid_ester          = 183;
  fg_phosphoric_acid_halide         = 184;
  fg_phosphoric_acid_amide          = 185;
  fg_thiophosphoric_acid_deriv      = 186;
  fg_thiophosphoric_acid            = 187;
  fg_thiophosphoric_acid_ester      = 188;
  fg_thiophosphoric_acid_halide     = 189;
  fg_thiophosphoric_acid_amide      = 190;
  fg_phosphonic_acid_deriv          = 191;
  fg_phosphonic_acid                = 192;
  fg_phosphonic_acid_ester          = 193;
  fg_phosphine                      = 194;
  fg_phosphinoxide                  = 195;
  fg_boronic_acid_deriv             = 196;
  fg_boronic_acid                   = 197;
  fg_boronic_acid_ester             = 198;
  fg_alkene                         = 199;
  fg_alkyne                         = 200;
  fg_aromatic                       = 201;
  fg_heterocycle                    = 202;
  fg_alpha_aminoacid                = 203;
  fg_alpha_hydroxyacid              = 204;

  max_atomicnum           = 118; // v0.4
  max_fragpath_length     = 8;   // v0.4
  min_fragpath_length     = 3;   // v0.4
  max_ringfragpath_length = 12;  // should be identical with ssr_vringsize
  max_fragstr_length      = 64;  // v0.4; must be at least 4 * max_ringfragpath_length
  hfpsize                 = 512; // v0.4
  n_hfpbits               = 2;   // v0.4; currently supported: 2 or 3
  max_n_matches           = 64;  // v0.4b; maximum number of matches of one and the same needle+haystack pair
  lang_en                 = 1;   // v0.5
  lang_de                 = 2;   // v0.5
  max_fgpos               = 127;  // v0.5 max. number of located funtional groups per fg type

type
  str2 = string[2];
  str3 = string[3];
  str4 = string[4];
  str5 = string[5];
  str8 = string[8];
  atom_rec  = record
                element : str2;
                atype : str3;
                x : single;
                y : single;
                z : single;
                formal_charge : integer;
                real_charge : single;
                Hexp : smallint;  // explicit H count
                Htot : smallint;  // total H count
                neighbor_count : integer;
                ring_count : integer;
                arom : boolean;
                q_arom : boolean;  // v0.3p potentially aromatic in a query structure
                stereo_care : boolean;   // new in v0.3d
	        heavy : boolean;         // new in v0.3l
		metal : boolean;         // new in v0.3l
		nvalences : integer;     // new in v0.3m
		tag : boolean;           // v0.3o
                nucleon_number : integer; // v0.3p
                radical_type : integer;   // v0.3p
              end;
  bond_rec  = record
                a1 : integer;
                a2 : integer;
                btype : char;
                ring_count : integer;
                arom : boolean;
                q_arom : boolean;  // v0.3p potentially aromatic in a query structure
                topo : shortint;   // new in v0.3d, see MDL file description
                stereo : shortint; // new in v0.3d
                mdl_stereo : shortint; // new in v0.3n
              end;
  ringpath_type = array[1..max_ringsize] of integer;
  matchpath_type = array[1..max_matchpath_length] of integer;

  atomlist = array[1..max_atoms] of atom_rec;
  bondlist = array[1..max_bonds] of bond_rec;
  ringlist = array[1..max_rings] of ringpath_type;
  neighbor_rec = array[1..max_neighbors] of integer;   // new in v0.2h
  fglist   = array[1..max_fg] of boolean;

  molbuftype = array[1..(max_atoms+max_bonds+8192)] of string;

  matchmatrix = array[1..max_neighbors,1..max_neighbors] of boolean;  // new in v0.2i
  global_matchmatrix = array[1..max_ndl_gmmsize,1..max_hst_gmmsize] of boolean;  // v0.4a

  molstat_rec = record
                  n_QA, n_QB, n_chg : integer;                  // number of query atoms, query bonds, charges
                  n_C1, n_C2, n_C : integer;                    // number of sp, sp2 hybridized, and total no. of carbons
                  n_CHB1p, n_CHB2p, n_CHB3p, n_CHB4 : integer;  // number of C atoms with at least 1, 2, 3 hetero bonds
                  n_O2, n_O3 : integer;                         // number of sp2 and sp3 oxygens
                  n_N1, n_N2, n_N3 : integer;                   // number of sp, sp2, and sp3 nitrogens
                  n_S, n_SeTe : integer;                        // number of sulfur atoms and selenium or tellurium atoms
                  n_F, n_Cl, n_Br, n_I : integer;               // number of fluorine, chlorine, bromine, iodine atoms
                  n_P, n_B : integer;                           // number of phosphorus and boron atoms
                  n_Met, n_X : integer;                         // number of metal and "other" atoms (not listed elsewhere); v0.3l
                  n_b1, n_b2, n_b3, n_bar : integer;            // number single, double, triple, and aromatic bonds
                  n_C1O, n_C2O, n_CN, n_XY : integer;           // number of C-O single bonds, C=O double bonds, CN bonds (any type), hetero/hetero bonds
                  n_r3, n_r4, n_r5, n_r6, n_r7, n_r8 : integer; // number of 3-, 4-, 5-, 6-, 7-, and 8-membered rings
                  n_r9, n_r10, n_r11, n_r12, n_r13p : integer;  // number of 9-, 10-, 11-, 12-, and 13plus-membered rings
                  n_rN, n_rN1, n_rN2, n_rN3p : integer;         // number of rings containing N (any number), 1 N, 2 N, and 3 N or more
                  n_rO, n_rO1, n_rO2p : integer;                // number of rings containing O (any number), 1 O, and 2 O or more
                  n_rS, n_rX, n_rAr, n_rBz  : integer;          // number of rings containing S (any number), any heteroatom (any number), 
                                                                // number of aromatic rings, number of benzene rings
                  n_br2p : integer;                             // number of bonds belonging to more than one ring (v0.3n)
                  {$IFDEF extended_molstat}
                  n_psg01, n_psg02, n_psg13, n_psg14 : integer; // number of atoms belonging to elements of group 1, 2, etc. 
                  n_psg15, n_psg16, n_psg17, n_psg18 : integer; // number of atoms belonging to elements of group 15, 16, etc. 
                  n_pstm, n_psla : integer;                     // number of transition metals, lanthanides/actinides
                  {$ENDIF}
                end;
  ringprop_rec  = record     // new in v0.3
                    size     : integer;
                    arom     : boolean;
                    envelope : boolean;
                  end;
  ringprop_type = array[1..max_rings] of ringprop_rec;
  p_3d          = record     // new in v0.3d
                    x : double;
                    y : double;
                    z : double;
                  end;
  chirpath_type = array[1..4] of integer;  // new in v0.3f
  connval_rec   = record                   // new in v0.3j
                    def : integer;         // better as longint for large molecules?
                    tmp : integer;
                  end;
  connval_type  = array[1..max_atoms] of connval_rec; // new in v0.3j 

  // new in v0.3q: periodic table lookup
  pt_rec = record
             el    : str2;
             am    : single;    // average isotopic mass; v0.4d  
             count : integer;   // v0.4d
           end;
  fragstr = string[max_fragstr_length];
  fragpath_type = array[1..max_ringfragpath_length] of integer;
  frag_rec = record
               size : integer;
               a_atomicnum : array[1..max_ringfragpath_length] of integer;
               b_code : array[1..max_ringfragpath_length] of char;
               ring : boolean;
             end;
  fgloctype = array[1..used_fg,0..max_fgpos] of smallint;  // v0.5
  
var
  progmode      : integer;
  progname      : string;
  i             : integer;  // general purpose index
  li            : longint;
  opt_none      : boolean;
  opt_verbose   : boolean;
  opt_text      : boolean;
  opt_text_de   : boolean;
  opt_code      : boolean;
  opt_bin       : boolean;
  opt_bitstring : boolean;
  opt_stdin     : boolean;
  opt_exact     : boolean;
  opt_debug     : boolean;
  opt_molout    : boolean;
  opt_molstat   : boolean;
  opt_molstat_X : boolean;
  opt_molstat_v : boolean;  // new in v0.4d
  opt_xmdlout   : boolean;
  opt_strict    : boolean;  // new in v0.2f
  opt_nochirality : boolean; // new in v0.4e (thanks to Robert Bruccoleri for this suggestion)
  opt_metalrings: boolean;  // new in v0.3
  opt_geom      : boolean;  // new in v0.3d
  opt_chiral    : boolean;  // new in v0.3f
  opt_iso       : boolean;  // new in v0.3p
  opt_chg       : boolean;  // new in v0.3p 
  opt_rad       : boolean;  // new in v0.3p
  opt_rs        : integer;  // new in v0.3i
  opt_fp        : boolean;  // new in v0.3m
  fpformat      : integer;  // new in v0.3m
  opt_hfp       : boolean;  // new in v0.4
  opt_matchnum  : boolean;  // new in v0.4a
  opt_matchnum1 : boolean;  // new in v0.4a
  opt_pos       : boolean;  // v0.5
  filetype : string;
  //molfile : text;
  molfilename : string;
  ndl_molfilename : string;
  molname : string;
  ndl_molname : string;
  tmp_molname : string;  // v0.3m
  molcomment  : string;
  n_atoms : integer;
  n_bonds : integer;
  n_rings : integer;    // the number of rings we determined ourselves
  n_countablerings : integer; // v0.3n; number of rings which are not envelope rings
  n_cmrings : integer;  // the number of rings we read from a (CheckMol-tweaked) MDL molfile
  n_charges : integer;  // number of charges    // unnecessary!!!!
  n_heavyatoms : integer;
  n_trueheavyatoms : integer;  // v0.4b
  n_heavybonds : integer;
  ndl_n_atoms : integer;
  ndl_n_bonds : integer;
  ndl_n_rings : integer;
  ndl_n_heavyatoms : integer;
  ndl_n_trueheavyatoms : integer;  // v0.4b
  ndl_n_heavybonds : integer;
  //cm_mdlmolfile  : boolean;
  found_arominfo : boolean;
  found_querymol : boolean;
  ndl_querymol : boolean;  // v0.3p
  tmp_n_atoms : integer;      // v0.3m
  tmp_n_bonds : integer;      // v0.3m
  tmp_n_rings : integer;      // v0.3m
  tmp_n_heavyatoms : integer; // v0.3m
  tmp_n_trueheavyatoms : integer;  // v0.4b
  tmp_n_heavybonds : integer; // v0.3m

  atom : ^atomlist;
  bond : ^bondlist;
  ring : ^ringlist;
  ringprop : ^ringprop_type;  // new in v0.3

  ndl_atom : ^atomlist;
  ndl_bond : ^bondlist;
  ndl_ring : ^ringlist;
  ndl_ringprop : ^ringprop_type;  // new in v0.3
  ndl_ref_atom : integer;  // since v0.3j as a global variable
  tmp_atom : ^atomlist;           // v0.3m
  tmp_bond : ^bondlist;           // v0.3m
  tmp_ring : ^ringlist;           // v0.3m
  tmp_ringprop : ^ringprop_type;  // v0.3m

  matchresult   : boolean;
  matchsummary  : boolean;        // v0.3o
  ndl_matchpath : matchpath_type;
  hst_matchpath : matchpath_type;

  fg   : fglist;
  atomtype : str4;
  newatomtype : str3;

  molbuf : ^molbuftype;
  molbufindex : integer;

  mol_in_queue : boolean;
  mol_count    : longint;

  molstat      : molstat_rec;
  ndl_molstat  : molstat_rec;
  tmp_molstat  : molstat_rec;  // v0.3m
  
  ringsearch_mode : integer;
  max_vringsize   : integer;  // for SSR ring search

  rfile : text;
  rfile_is_open : boolean;
  mol_OK        : boolean;  // new in v0.2i

  n_ar          : integer;  // new in v0.3
  prev_n_ar     : integer;  // new in v0.3
  ez_search     : boolean;  // new in v0.3d
  rs_search     : boolean;  // new in v0.3f
  ez_flag       : boolean;  // new in v0.3f
  chir_flag     : boolean;  // new in v0.3f
  rs_strict     : boolean;  // new in v0.3j

  n_Ctot, n_Otot, n_Ntot             : integer;  // new in v0.3g
  ndl_n_Ctot, ndl_n_Otot, ndl_n_Ntot : integer;  // new in v0.3g
  tmp_n_Ctot, tmp_n_Otot, tmp_n_Ntot : integer;  // new in v0.3m
  ether_generic : boolean;   // v0.3j
  amine_generic : boolean;   // v0.3j  
  hydroxy_generic : boolean; // v0.3j
  cv            : ^connval_type;  // new in v0.3j
  fpdecimal     : int64;  // v0.3m
  fpincrement   : int64;
  fpindex       : integer;
  fp_exacthit   : boolean;
  fp_exactblock : boolean;
  tmfcode       : integer;  // v0.3m; version number for tweaked MDL molfiles (tweaklevel)
  tmfmismatch   : boolean;  // v0.3m; rely on tweaked MDL molfiles only if same level
  auto_ssr      : boolean;  // v0.3n; indicates that SAR -> SSR fallback has happened
  recursion_level : longint; // v0.3p
  recursion_depth : longint; // v0.3p
  keep_DT         : boolean; // v0.3p  molfile output: Deuterium/Tritium as D/T or 2H/3H

  // new in v0.4: periodic table lookup, fragments for hashed fingerprints
  pt : array[1..max_atomicnum] of pt_rec;
  hfp : array[1..hfpsize] of boolean;
  hfpformat      : integer;  // end of new v0.4 variables
  // new in v0.4a: global match matrix (+ cumulative copy), various auxiliary variables
  gmm            : ^global_matchmatrix;
  gmm_total      : ^global_matchmatrix;
  use_gmm        : boolean;
  valid_gmm      : boolean;
  overall_match  : boolean;
  tmp_tag        : array[1..max_atoms] of boolean;
  found_untagged : boolean;  // end of new v0.4a variables
  n_matches      : integer;  // v0.4b
  match_string   : ansistring;  // v0.4b
  sep_label      : string;  // v0.4b
  mol_formula    : string;  // v0.4d
  mol_weight     : single;  // v0.4d
  fglang         : integer;  // v0.5
  fgloc          : ^fgloctype;
  
//============================= auxiliary functions & procedures


function mkfglabel(fgnum:integer;lang:integer):string;  // v0.5
var
  f : string;
begin
  if ((fgnum < 0) or (fgnum > used_fg) or (lang < 0)) then
    begin
      mkfglabel := '';
      exit;
    end;
  f := '';
  if (lang = 0) then   // lang 0: 8-character fg codes
    begin
      case fgnum of
        1   :  f := '000000T2';     // cation                                                    
        2   :  f := '000000T1';     // anion                                                     
        3   :  f := 'C2O10000';     // carbonyl compound                                         
        4   :  f := 'C2O1H000';     // aldehyde                                                  
        5   :  f := 'C2O1C000';     // ketone                                                    
        6   :  f := 'C2S10000';     // thiocarbonyl compound                                     
        7   :  f := 'C2S1H000';     // thioaldehyde                                              
        8   :  f := 'C2S1C000';     // thioketone                                                
        9   :  f := 'C2N10000';     // imine                                                     
        10  :  f := 'C2N1N000';     // hydrazone                                                 
        11  :  f := 'C2NNC4ON';     // semicarbazone                                             
        12  :  f := 'C2NNC4SN';     // thiosemicarbazone                                         
        13  :  f := 'C2N1OH00';     // oxime                                                     
        14  :  f := 'C2N1OC00';     // oxime ether                                               
        15  :  f := 'C3OC0000';     // ketene                                                    
        16  :  f := 'C3OCC000';     // ketene acetal or derivative                               
        17  :  f := 'C2O2H200';     // carbonyl hydrate                                          
        18  :  f := 'C2O2HC00';     // hemiacetal                                                
        19  :  f := 'C2O2CC00';     // acetal                                                    
        20  :  f := 'C2NOHC10';     // hemiaminal                                                
        21  :  f := 'C2N2CC10';     // aminal                                                    
        22  :  f := 'C2NSHC10';     // hemithioaminal                                            
        23  :  f := 'C2S2CC00';     // thioacetal                                                
        24  :  f := 'C2CNH000';     // enamine                                                   
        25  :  f := 'C2COH000';     // enol                                                      
        26  :  f := 'C2COH000';     // enol ether                                                
        27  :  f := 'O1H00000';     // hydroxy compound                                          
        28  :  f := 'O1H0C000';     // alcohol                                                   
        29  :  f := 'O1H1C000';     // primary alcohol                                           
        30  :  f := 'O1H2C000';     // secondary alcohol                                         
        31  :  f := 'O1H3C000';     // tertiary alcohol                                          
        32  :  f := 'O1H0CO1H';     // 1,2-diol                                                  
        33  :  f := 'O1H0CN1C';     // 1,2-aminoalcohol                                          
        34  :  f := 'O1H1A000';     // phenol or hydroxyhetarene                                 
        35  :  f := 'O1H2A000';     // 1,2-diphenol                                              
        36  :  f := 'C2COH200';     // enediol                                                   
        37  :  f := 'O1C00000';     // ether                                                     
        38  :  f := 'O1C0CC00';     // dialkyl ether                                             
        39  :  f := 'O1C0CA00';     // alkyl aryl ether                                          
        40  :  f := 'O1C0AA00';     // diaryl ether                                              
        41  :  f := 'S1C00000';     // thioether                                                 
        42  :  f := 'S1S1C000';     // disulfide                                                 
        43  :  f := 'O1O1C000';     // peroxide                                                  
        44  :  f := 'O1O1H000';     // hydroperoxide                                             
        45  :  f := 'N1N10000';     // hydrazine derivative                                      
        46  :  f := 'N1O1H000';     // hydroxylamine                                             
        47  :  f := 'N1C00000';     // amine                                                     
        48  :  f := 'N1C10000';     // primary amine                                             
        49  :  f := 'N1C1C000';     // primary aliphatic amine (alkylamine)                      
        50  :  f := 'N1C1A000';     // primary aromatic amine                                    
        51  :  f := 'N1C20000';     // secondary amine                                           
        52  :  f := 'N1C2CC00';     // secondary aliphatic amine (dialkylamine)                  
        53  :  f := 'N1C2AC00';     // secondary aliphatic/aromatic amine (alkylarylamine)       
        54  :  f := 'N1C2AA00';     // secondary aromatic amine (diarylamine)                    
        55  :  f := 'N1C30000';     // tertiary amine                                            
        56  :  f := 'N1C3CC00';     // tertiary aliphatic amine (trialkylamine)                  
        57  :  f := 'N1C3AC00';     // tertiary aliphatic/aromatic amine (alkylarylamine)        
        58  :  f := 'N1C3AA00';     // tertiary aromatic amine (triarylamine)                    
        59  :  f := 'N1C400T2';     // quaternary ammonium salt                                  
        60  :  f := 'N0O10000';     // N-oxide                                                   
        61  :  f := 'XX000000';     // halogen derivative                                        
        62  :  f := 'XX00C000';     // alkyl halide                                              
        63  :  f := 'XF00C000';     // alkyl fluoride                                            
        64  :  f := 'XC00C000';     // alkyl chloride                                            
        65  :  f := 'XB00C000';     // alkyl bromide                                             
        66  :  f := 'XI00C000';     // alkyl iodide                                              
        67  :  f := 'XX00A000';     // aryl halide                                               
        68  :  f := 'XF00A000';     // aryl fluoride                                             
        69  :  f := 'XC00A000';     // aryl chloride                                             
        70  :  f := 'XB00A000';     // aryl bromide                                              
        71  :  f := 'XI00A000';     // aryl iodide                                               
        72  :  f := '000000MX';     // organometallic compound                                   
        73  :  f := '000000ML';     // organolithium compound                                    
        74  :  f := '000000MM';     // organomagnesium compound                                  
        75  :  f := 'C3O20000';     // carboxylic acid derivative                                
        76  :  f := 'C3O2H000';     // carboxylic acid                                           
        77  :  f := 'C3O200T1';     // carboxylic acid salt                                      
        78  :  f := 'C3O2C000';     // carboxylic acid ester                                     
        79  :  f := 'C3O2CZ00';     // lactone                                                      
        80  :  f := 'C3ONC000';     // carboxylic acid amide                                     
        81  :  f := 'C3ONC100';     // primary carboxylic acid amide                             
        82  :  f := 'C3ONC200';     // secondary carboxylic acid amide                           
        83  :  f := 'C3ONC300';     // tertiary carboxylic acid amide                            
        84  :  f := 'C3ONCZ00';     // lactam                                                        
        85  :  f := 'C3ONN100';     // carboxylic acid hydrazide                                 
        86  :  f := 'C3ONN200';     // carboxylic acid azide                                     
        87  :  f := 'C3ONOH00';     // hydroxamic acid                                           
        88  :  f := 'C3N2H000';     // carboxylic acid amidine                                   
        89  :  f := 'C3NNN100';     // carboxylic acid amidrazone                                
        90  :  f := 'C3N00000';     // carbonitrile                                              
        91  :  f := 'C3OXX000';     // acyl halide                                               
        92  :  f := 'C3OXF000';     // acyl fluoride                                             
        93  :  f := 'C3OXC000';     // acyl chloride                                             
        94  :  f := 'C3OXB000';     // acyl bromide                                              
        95  :  f := 'C3OXI000';     // acyl iodide                                               
        96  :  f := 'C2OC3N00';     // acyl cyanide                                              
        97  :  f := 'C3NOC000';     // imido ester                                               
        98  :  f := 'C3NXX000';     // imidoyl halide                                            
        99  :  f := 'C3SO0000';     // thiocarboxylic acid derivative                            
        100 :  f := 'C3SOH000';     // thiocarboxylic acid                                       
        101 :  f := 'C3SOC000';     // thiocarboxylic acid ester                                 
        102 :  f := 'C3SOCZ00';     // thiolactone                                               
        103 :  f := 'C3SNH000';     // thiocarboxylic acid amide                                 
        104 :  f := 'C3SNCZ00';     // thiolactam                                                
        105 :  f := 'C3NSC000';     // imidothioester                                            
        106 :  f := 'C3ONAZ00';     // oxo(het)arene                                             
        107 :  f := 'C3SNAZ00';     // thioxo(het)arene                                          
        108 :  f := 'C3NNAZ00';     // imino(het)arene                                           
        109 :  f := 'C3O30000';     // orthocarboxylic acid derivative                           
        110 :  f := 'C3O3C000';     // orthoester                                                
        111 :  f := 'C3O3NC00';     // amide acetal                                              
        112 :  f := 'C3O2C3O2';     // carboxylic acid anhydride                                 
        113 :  f := 'C3ONC000';     // carboxylic acid imide                                     
        114 :  f := 'C3ONCH10';     // carboxylic acid imide, N-unsubstituted                    
        115 :  f := 'C3ONCC10';     // carboxylic acid imide, N-substituted                      
        116 :  f := 'C4000000';     // CO2 derivative (general)                                  
        117 :  f := 'C4O30000';     // carbonic acid derivative                                  
        118 :  f := 'C4O3C100';     // carbonic acid monoester                                   
        119 :  f := 'C4O3C200';     // carbonic acid diester                                     
        120 :  f := 'C4O3CX00';     // carbonic acid ester halide (alkyl/aryl haloformate)       
        121 :  f := 'C4SO0000';     // thiocarbonic acid derivative                              
        122 :  f := 'C4SOC100';     // thiocarbonic acid monoester                               
        123 :  f := 'C4SOC200';     // thiocarbonic acid diester                                 
        124 :  f := 'C4SOX_00';     // thiocarbonic acid ester halide (alkyl/aryl halothioformate
        125 :  f := 'C4O2N000';     // carbamic acid derivative                                  
        126 :  f := 'C4O2NH00';     // carbamic acid                                             
        127 :  f := 'C4O2NC00';     // carbamic acid ester (urethane)                            
        128 :  f := 'C4O2NX00';     // carbamic acid halide (haloformic acid amide)              
        129 :  f := 'C4SN0000';     // thiocarbamic acid derivative                              
        130 :  f := 'C4SNOH00';     // thiocarbamic acid                                         
        131 :  f := 'C4SNOC00';     // thiocarbamic acid ester                                   
        132 :  f := 'C4SNXX00';     // thiocarbamic acid halide (halothioformic acid amide)      
        133 :  f := 'C4O1N200';     // urea                                                           
        134 :  f := 'C4N2O100';     // isourea                                                     
        135 :  f := 'C4S1N200';     // thiourea                                                   
        136 :  f := 'C4N2S100';     // isothiourea                                               
        137 :  f := 'C4N30000';     // guanidine                                                 
        138 :  f := 'C4ON2N00';     // semicarbazide                                             
        139 :  f := 'C4SN2N00';     // thiosemicarbazide                                         
        140 :  f := 'N4N20000';     // azide                                                         
        141 :  f := 'N2N10000';     // azo compound                                              
        142 :  f := 'N3N100T2';     // diazonium salt                                            
        143 :  f := 'N3C10000';     // isonitrile                                                
        144 :  f := 'C4NO1000';     // cyanate                                                     
        145 :  f := 'C4NO2000';     // isocyanate                                                
        146 :  f := 'C4NS1000';     // thiocyanate                                               
        147 :  f := 'C4NS2000';     // isothiocyanate                                            
        148 :  f := 'C4N20000';     // carbodiimide                                              
        149 :  f := 'N2O10000';     // nitroso compound                                          
        150 :  f := 'N4O20000';     // nitro compound                                            
        151 :  f := 'N3O20000';     // nitrite                                                     
        152 :  f := 'N4O30000';     // nitrate                                                     
        153 :  f := 'S6O00000';     // sulfuric acid derivative                                  
        154 :  f := 'S6O4H000';     // sulfuric acid                                             
        155 :  f := 'S6O4HC00';     // sulfuric acid monoester                                   
        156 :  f := 'S6O4CC00';     // sulfuric acid diester                                     
        157 :  f := 'S6O3NC00';     // sulfuric acid amide ester                                 
        158 :  f := 'S6O3N100';     // sulfuric acid amide                                       
        159 :  f := 'S6O2N200';     // sulfuric acid diamide                                     
        160 :  f := 'S6O3XX00';     // sulfuryl halide                                           
        161 :  f := 'S5O00000';     // sulfonic acid derivative                                  
        162 :  f := 'S5O3H000';     // sulfonic acid                                             
        163 :  f := 'S5O3C000';     // sulfonic acid ester                                       
        164 :  f := 'S5O2N000';     // sulfonamide                                               
        165 :  f := 'S5O2XX00';     // sulfonyl halide                                           
        166 :  f := 'S4O20000';     // sulfone                                                     
        167 :  f := 'S2O10000';     // sulfoxide                                                 
        168 :  f := 'S3O00000';     // sulfinic acid derivative                                  
        169 :  f := 'S3O2H000';     // sulfinic acid                                             
        170 :  f := 'S3O2C000';     // sulfinic acid ester                                       
        171 :  f := 'S3O1XX00';     // sulfinic acid halide                                      
        172 :  f := 'S3O1N000';     // sulfinic acid amide                                       
        173 :  f := 'S1O00000';     // sulfenic acid derivative                                  
        174 :  f := 'S1O1H000';     // sulfenic acid                                             
        175 :  f := 'S1O1C000';     // sulfenic acid ester                                       
        176 :  f := 'S1O0XX00';     // sulfenic acid halide                                      
        177 :  f := 'S1O0N100';     // sulfenic acid amide                                       
        178 :  f := 'S1H10000';     // thiol (sulfanyl compound)                                 
        179 :  f := 'S1H1C000';     // alkylthiol                                                
        180 :  f := 'S1H1A000';     // arylthiol (thiophenol)                                    
        181 :  f := 'P5O0H000';     // phosphoric acid derivative                                
        182 :  f := 'P5O4H200';     // phosphoric acid                                           
        183 :  f := 'P5O4HC00';     // phosphoric acid ester                                     
        184 :  f := 'P5O3HX00';     // phosphoric acid halide                                    
        185 :  f := 'P5O3HN00';     // phosphoric acid amide                                     
        186 :  f := 'P5O0S000';     // thiophosphoric acid derivative                            
        187 :  f := 'P5O3SH00';     // thiophosphoric acid                                       
        188 :  f := 'P5O3SC00';     // thiophosphoric acid ester                                 
        189 :  f := 'P5O2SX00';     // thiophosphoric acid halide                                
        190 :  f := 'P5O2SN00';     // thiophosphoric acid amide                                 
        191 :  f := 'P4O30000';     // phosphonic acid derivative                                
        192 :  f := 'P4O3H000';     // phosphonic acid                                           
        193 :  f := 'P4O3C000';     // phosphonic acid ester                                     
        194 :  f := 'P3000000';     // phosphine                                                 
        195 :  f := 'P2O00000';     // phosphine oxide                                           
        196 :  f := 'B2O20000';     // boronic acid derivative                                   
        197 :  f := 'B2O2H000';     // boronic acid                                              
        198 :  f := 'B2O2C000';     // boronic acid ester                                        
        199 :  f := '000C2C00';     // alkene                                                       
        200 :  f := '000C3C00';     // alkyne                                                       
        201 :  f := '0000A000';     // aromatic compound                                         
        202 :  f := '0000CZ00';     // heterocyclic compound                                     
        203 :  f := 'C3O2HN1C';     // alpha-aminoacid                                           
        204 :  f := 'C3O2HO1H';     // alpha-hydroxyacid                                         
      end;  // case
    end;   // if lang = 0
  if (lang = lang_en) then   // lang 1: english
    begin
      case fgnum of
        1   :  f := 'cation';
        2   :  f := 'anion';
        3   :  f := 'carbonyl compound';
        4   :  f := 'aldehyde';
        5   :  f := 'ketone';
        6   :  f := 'thiocarbonyl compound';
        7   :  f := 'thioaldehyde';
        8   :  f := 'thioketone';
        9   :  f := 'imine';
        10  :  f := 'hydrazone';
        11  :  f := 'semicarbazone';
        12  :  f := 'thiosemicarbazone';
        13  :  f := 'oxime';
        14  :  f := 'oxime ether';
        15  :  f := 'ketene';
        16  :  f := 'ketene acetal or derivative';
        17  :  f := 'carbonyl hydrate';
        18  :  f := 'hemiacetal';
        19  :  f := 'acetal';
        20  :  f := 'hemiaminal';
        21  :  f := 'aminal';
        22  :  f := 'hemithioaminal';
        23  :  f := 'thioacetal';
        24  :  f := 'enamine';
        25  :  f := 'enol';
        26  :  f := 'enol ether';
        27  :  f := 'hydroxy compound';
        28  :  f := 'alcohol';
        29  :  f := 'primary alcohol';
        30  :  f := 'secondary alcohol';
        31  :  f := 'tertiary alcohol';
        32  :  f := '1,2-diol';
        33  :  f := '1,2-aminoalcohol';
        34  :  f := 'phenol or hydroxyhetarene';
        35  :  f := '1,2-diphenol';
        36  :  f := 'enediol';
        37  :  f := 'ether';
        38  :  f := 'dialkyl ether';
        39  :  f := 'alkyl aryl ether';
        40  :  f := 'diaryl ether';
        41  :  f := 'thioether';
        42  :  f := 'disulfide';
        43  :  f := 'peroxide';
        44  :  f := 'hydroperoxide';
        45  :  f := 'hydrazine derivative';
        46  :  f := 'hydroxylamine';
        47  :  f := 'amine';
        48  :  f := 'primary amine';
        49  :  f := 'primary aliphatic amine (alkylamine)';
        50  :  f := 'primary aromatic amine';
        51  :  f := 'secondary amine';
        52  :  f := 'secondary aliphatic amine (dialkylamine)';
        53  :  f := 'secondary aliphatic/aromatic amine (alkylarylamine)';
        54  :  f := 'secondary aromatic amine (diarylamine)';
        55  :  f := 'tertiary amine';
        56  :  f := 'tertiary aliphatic amine (trialkylamine)';
        57  :  f := 'tertiary aliphatic/aromatic amine (alkylarylamine)';
        58  :  f := 'tertiary aromatic amine (triarylamine)';
        59  :  f := 'quaternary ammonium salt';
        60  :  f := 'N-oxide';
        61  :  f := 'halogen derivative';
        62  :  f := 'alkyl halide';
        63  :  f := 'alkyl fluoride';
        64  :  f := 'alkyl chloride'; 
        65  :  f := 'alkyl bromide';
        66  :  f := 'alkyl iodide';
        67  :  f := 'aryl halide'; 
        68  :  f := 'aryl fluoride';
        69  :  f := 'aryl chloride';
        70  :  f := 'aryl bromide';
        71  :  f := 'aryl iodide'; 
        72  :  f := 'organometallic compound';
        73  :  f := 'organolithium compound';
        74  :  f := 'organomagnesium compound';
        75  :  f := 'carboxylic acid derivative';
        76  :  f := 'carboxylic acid';
        77  :  f := 'carboxylic acid salt';
        78  :  f := 'carboxylic acid ester';
        79  :  f := 'lactone';     
        80  :  f := 'carboxylic acid amide';
        81  :  f := 'primary carboxylic acid amide';
        82  :  f := 'secondary carboxylic acid amide';
        83  :  f := 'tertiary carboxylic acid amide';
        84  :  f := 'lactam';      
        85  :  f := 'carboxylic acid hydrazide';
        86  :  f := 'carboxylic acid azide';
        87  :  f := 'hydroxamic acid';
        88  :  f := 'carboxylic acid amidine';
        89  :  f := 'carboxylic acid amidrazone';
        90  :  f := 'carbonitrile';
        91  :  f := 'acyl halide'; 
        92  :  f := 'acyl fluoride';
        93  :  f := 'acyl chloride';
        94  :  f := 'acyl bromide';
        95  :  f := 'acyl iodide'; 
        96  :  f := 'acyl cyanide';
        97  :  f := 'imido ester'; 
        98  :  f := 'imidoyl halide';
        99  :  f := 'thiocarboxylic acid derivative';
        100 :  f := 'thiocarboxylic acid';
        101 :  f := 'thiocarboxylic acid ester';
        102 :  f := 'thiolactone';
        103 :  f := 'thiocarboxylic acid amide';
        104 :  f := 'thiolactam'; 
        105 :  f := 'imidothioester';
        106 :  f := 'oxo(het)arene';
        107 :  f := 'thioxo(het)arene';
        108 :  f := 'imino(het)arene';
        109 :  f := 'orthocarboxylic acid derivative';
        110 :  f := 'orthoester'; 
        111 :  f := 'amide acetal';
        112 :  f := 'carboxylic acid anhydride';
        113 :  f := 'carboxylic acid imide';
        114 :  f := 'carboxylic acid imide, N-unsubstituted';
        115 :  f := 'carboxylic acid imide, N-substituted';
        116 :  f := 'CO2 derivative (general)';
        117 :  f := 'carbonic acid derivative';
        118 :  f := 'carbonic acid monoester';
        119 :  f := 'carbonic acid diester';
        120 :  f := 'carbonic acid ester halide (alkyl/aryl haloformate)'; 
        121 :  f := 'thiocarbonic acid derivative';
        122 :  f := 'thiocarbonic acid monoester';
        123 :  f := 'thiocarbonic acid diester';
        124 :  f := 'thiocarbonic acid ester halide (alkyl/aryl halothioformate';
        125 :  f := 'carbamic acid derivative';
        126 :  f := 'carbamic acid';
        127 :  f := 'carbamic acid ester (urethane)';
        128 :  f := 'carbamic acid halide (haloformic acid amide)';
        129 :  f := 'thiocarbamic acid derivative';
        130 :  f := 'thiocarbamic acid';
        131 :  f := 'thiocarbamic acid ester';
        132 :  f := 'thiocarbamic acid halide (halothioformic acid amide)';
        133 :  f := 'urea';       
        134 :  f := 'isourea';    
        135 :  f := 'thiourea';   
        136 :  f := 'isothiourea';
        137 :  f := 'guanidine';  
        138 :  f := 'semicarbazide';
        139 :  f := 'thiosemicarbazide';
        140 :  f := 'azide';      
        141 :  f := 'azo compound';
        142 :  f := 'diazonium salt';
        143 :  f := 'isonitrile'; 
        144 :  f := 'cyanate';    
        145 :  f := 'isocyanate'; 
        146 :  f := 'thiocyanate';
        147 :  f := 'isothiocyanate';
        148 :  f := 'carbodiimide';
        149 :  f := 'nitroso compound';
        150 :  f := 'nitro compound';
        151 :  f := 'nitrite';    
        152 :  f := 'nitrate';    
        153 :  f := 'sulfuric acid derivative';
        154 :  f := 'sulfuric acid';
        155 :  f := 'sulfuric acid monoester';
        156 :  f := 'sulfuric acid diester';
        157 :  f := 'sulfuric acid amide ester';
        158 :  f := 'sulfuric acid amide';
        159 :  f := 'sulfuric acid diamide';
        160 :  f := 'sulfuryl halide';
        161 :  f := 'sulfonic acid derivative';
        162 :  f := 'sulfonic acid';
        163 :  f := 'sulfonic acid ester';
        164 :  f := 'sulfonamide';
        165 :  f := 'sulfonyl halide';
        166 :  f := 'sulfone';    
        167 :  f := 'sulfoxide';  
        168 :  f := 'sulfinic acid derivative';
        169 :  f := 'sulfinic acid';
        170 :  f := 'sulfinic acid ester';
        171 :  f := 'sulfinic acid halide';
        172 :  f := 'sulfinic acid amide';
        173 :  f := 'sulfenic acid derivative';
        174 :  f := 'sulfenic acid';
        175 :  f := 'sulfenic acid ester';
        176 :  f := 'sulfenic acid halide';
        177 :  f := 'sulfenic acid amide';
        178 :  f := 'thiol (sulfanyl compound)';
        179 :  f := 'alkylthiol'; 
        180 :  f := 'arylthiol (thiophenol)';
        181 :  f := 'phosphoric acid derivative';
        182 :  f := 'phosphoric acid';
        183 :  f := 'phosphoric acid ester';
        184 :  f := 'phosphoric acid halide';
        185 :  f := 'phosphoric acid amide';
        186 :  f := 'thiophosphoric acid derivative';
        187 :  f := 'thiophosphoric acid';
        188 :  f := 'thiophosphoric acid ester';
        189 :  f := 'thiophosphoric acid halide';
        190 :  f := 'thiophosphoric acid amide';
        191 :  f := 'phosphonic acid derivative';
        192 :  f := 'phosphonic acid';
        193 :  f := 'phosphonic acid ester';
        194 :  f := 'phosphine';  
        195 :  f := 'phosphine oxide';
        196 :  f := 'boronic acid derivative';
        197 :  f := 'boronic acid';
        198 :  f := 'boronic acid ester';
        199 :  f := 'alkene';     
        200 :  f := 'alkyne';     
        201 :  f := 'aromatic compound';
        202 :  f := 'heterocyclic compound';
        203 :  f := 'alpha-aminoacid';
        204 :  f := 'alpha-hydroxyacid';
      end;  // case
    end;   // if lang = lang_en
  if (lang = lang_de) then   // lang 2: german
    begin
      case fgnum of
        1   :  f := 'Kation';
        2   :  f := 'Anion';
        3   :  f := 'Carbonylverbindung';
        4   :  f := 'Aldehyd';
        5   :  f := 'Keton';
        6   :  f := 'Thiocarbonylverbindung';
        7   :  f := 'Thioaldehyd';
        8   :  f := 'Thioketon';
        9   :  f := 'Imin';
        10  :  f := 'Hydrazon';
        11  :  f := 'Semicarbazon';
        12  :  f := 'Thiosemicarbazon';
        13  :  f := 'Oxim';
        14  :  f := 'Oximether';
        15  :  f := 'Keten';
        16  :  f := 'Keten-Acetal oder Derivat';
        17  :  f := 'Carbonyl-Hydrat';
        18  :  f := 'Halbacetal';
        19  :  f := 'Halbacetal';
        20  :  f := 'Halbaminal';
        21  :  f := 'Aminal';
        22  :  f := 'Thiohalbaminal';
        23  :  f := 'Thioacetal';
        24  :  f := 'Enamin';
        25  :  f := 'Enol';
        26  :  f := 'Enolether';
        27  :  f := 'Hydroxy-Verbindung';
        28  :  f := 'Alkohol';
        29  :  f := 'primärer Alkohol';
        30  :  f := 'sekundärer Alkohol';
        31  :  f := 'tertiärer Alkohol';
        32  :  f := '1,2-Diol';
        33  :  f := '1,2-Aminoalkohol';
        34  :  f := 'Phenol oder Hydroxyhetaren';
        35  :  f := '1,2-Diphenol';
        36  :  f := 'Endiol';
        37  :  f := 'Ether';
        38  :  f := 'Dialkylether';
        39  :  f := 'Alkylarylether';
        40  :  f := 'Diarylether';
        41  :  f := 'Thioether';
        42  :  f := 'Disulfid';
        43  :  f := 'Peroxid';
        44  :  f := 'Hydroperoxid';
        45  :  f := 'Hydrazin-Derivat';
        46  :  f := 'Hydroxylamin';
        47  :  f := 'Amin';
        48  :  f := 'primäres Amin';
        49  :  f := 'primäres aliphatisches Amin (Alkylamin)';
        50  :  f := 'primäres aromatisches Amin';
        51  :  f := 'sekundäres Amin';
        52  :  f := 'sekundäres aliphatisches Amin (Dialkylamin)';
        53  :  f := 'sekundäres aliphatisches Amin (Dialkylamin)';
        54  :  f := 'sekundäres aromatisches Amin (Diarylamin)';
        55  :  f := 'tertiäres Amin';
        56  :  f := 'tertiäres aliphatisches Amin (Trialkylamin)';
        57  :  f := 'tertiäres aliphatisches/aromatisches Amin (Alkylarylamin)';
        58  :  f := 'tertiäres aromatisches Amin (Triarylamin)';
        59  :  f := 'quartäres Ammoniumsalz';
        60  :  f := 'N-Oxid';
        61  :  f := 'Halogenverbindung';
        62  :  f := 'Alkylhalogenid';
        63  :  f := 'Alkylfluorid';
        64  :  f := 'Alkylchlorid';
        65  :  f := 'Alkylbromid';
        66  :  f := 'Alkyliodid';
        67  :  f := 'Arylhalogenid';
        68  :  f := 'Arylfluorid';
        69  :  f := 'Arylchlorid';
        70  :  f := 'Arylbromid';
        71  :  f := 'Aryliodid';
        72  :  f := 'Organometall-Verbindung';
        73  :  f := 'Organolithium-Verbindung';
        74  :  f := 'Organomagnesium-Verbindung';
        75  :  f := 'Carbonsäure-Derivat';
        76  :  f := 'Carbonsäure';
        77  :  f := 'Carbonsäuresalz';
        78  :  f := 'Carbonsäureester';
        79  :  f := 'Lacton';
        80  :  f := 'Carbonsäureamid';
        81  :  f := 'primäres Carbonsäureamid';
        82  :  f := 'sekundäres Carbonsäureamid';
        83  :  f := 'tertiäres Carbonsäureamid';
        84  :  f := 'Lactam';
        85  :  f := 'Carbonsäurehydrazid';
        86  :  f := 'Carbonsäureazid';
        87  :  f := 'Hydroxamsäure';
        88  :  f := 'Carbonsäureamidin';
        89  :  f := 'Carbonsäureamidrazon';
        90  :  f := 'Carbonitril';
        91  :  f := 'Acylhalogenid';
        92  :  f := 'Acylfluorid';
        93  :  f := 'Acylchlorid';
        94  :  f := 'Acylbromid';
        95  :  f := 'Acyliodid';
        96  :  f := 'Acylcyanid';
        97  :  f := 'Imidoester';
        98  :  f := 'Imidoylhalogenid';
        99  :  f := 'Thiocarbonsäure-Derivat';
        100 :  f := 'Thiocarbonsäure';
        101 :  f := 'Thiocarbonsäureester';
        102 :  f := 'Thiolacton';
        103 :  f := 'Thiocarbonsäureamid';
        104 :  f := 'Thiolactam';
        105 :  f := 'Imidothioester';
        106 :  f := 'Oxo(het)aren';
        107 :  f := 'Thioxo(het)aren';
        108 :  f := 'Imino(het)aren';
        109 :  f := 'Orthocarbonsäure-Derivat';
        110 :  f := 'Orthoester';
        111 :  f := 'Amidacetal';
        112 :  f := 'Carbonsäureanhydrid';
        113 :  f := 'Carbonsäureimid';
        114 :  f := 'Carbonsäureimid, N-unsubstituiert';
        115 :  f := 'Carbonsäureimid, N-substituiert';
        116 :  f := 'CO2-Derivat (allgemein)';
        117 :  f := 'Kohlensäure-Derivat';
        118 :  f := 'Kohlensäuremonoester';
        119 :  f := 'Kohlensäurediester';
        120 :  f := 'Kohlensäureesterhalogenid (Alkyl/Aryl-Halogenformiat)';
        121 :  f := 'Thiokohlensäure-Derivat';
        122 :  f := 'Thiokohlensäuremonoester';
        123 :  f := 'Thiokohlensäurediester';
        124 :  f := 'Thiokohlensäureesterhalogenid (Alkyl/Aryl-Halogenthioformiat)';
        125 :  f := 'Carbaminsäure-Derivat';
        126 :  f := 'Carbaminsäure';
        127 :  f := 'Carbaminsäureester (Urethan)';
        128 :  f := 'Carbaminsäurehalogenid (Halogenformamid)';
        129 :  f := 'Thiocarbaminsäure-Derivat';
        130 :  f := 'Thiocarbaminsäure';
        131 :  f := 'Thiocarbaminsäureester';
        132 :  f := 'Thiocarbaminsäurehalogenid (Halogenthioformamid)';
        133 :  f := 'Harnstoff';
        134 :  f := 'Isoharnstoff';
        135 :  f := 'Thioharnstoff';
        136 :  f := 'Isothioharnstoff';
        137 :  f := 'Guanidin';
        138 :  f := 'Semicarbazid';
        139 :  f := 'Thiosemicarbazid';
        140 :  f := 'Azid';
        141 :  f := 'Azoverbindung';
        142 :  f := 'Diazoniumsalz';
        143 :  f := 'Isonitril';
        144 :  f := 'Cyanat';
        145 :  f := 'Isocyanat';
        146 :  f := 'Thiocyanat';
        147 :  f := 'Isothiocyanat';
        148 :  f := 'Carbodiimid';
        149 :  f := 'Nitroso-Verbindung';
        150 :  f := 'Nitro-Verbindung';
        151 :  f := 'Nitrit';
        152 :  f := 'Nitrat';
        153 :  f := 'Schwefelsäure-Derivat';
        154 :  f := 'Schwefelsäure';
        155 :  f := 'Schwefelsäuremonoester';
        156 :  f := 'Schwefelsäurediester';
        157 :  f := 'Schwefelsäureamidester';
        158 :  f := 'Schwefelsäureamid';
        159 :  f := 'Schwefelsäurediamid';
        160 :  f := 'Sulfurylhalogenid';
        161 :  f := 'Sulfonsäure-Derivat';
        162 :  f := 'Sulfonsäure';
        163 :  f := 'Sulfonsäureester';
        164 :  f := 'Sulfonamid';
        165 :  f := 'Sulfonylhalogenid';
        166 :  f := 'Sulfon';
        167 :  f := 'Sulfoxid';
        168 :  f := 'Sulfinsäure-Derivat';
        169 :  f := 'Sulfinsäure';
        170 :  f := 'Sulfinsäureester';
        171 :  f := 'Sulfinsäurehalogenid';
        172 :  f := 'Sulfinsäureamid';
        173 :  f := 'Sulfensäure-Derivat';
        174 :  f := 'Sulfensäure';
        175 :  f := 'Sulfensäureester';
        176 :  f := 'Sulfensäurehalogenid';
        177 :  f := 'Sulfensäureamid';
        178 :  f := 'Thiol (Sulfanyl-Verbindung, Mercaptan)';
        179 :  f := 'Alkylthiol';
        180 :  f := 'Arylthiol (Thiophenol)';
        181 :  f := 'Phosphorsäure-Derivat';
        182 :  f := 'Phosphorsäure';
        183 :  f := 'Phosphorsäureester';
        184 :  f := 'Phosphorsäurehalogenid';
        185 :  f := 'Phosphorsäureamid';
        186 :  f := 'Thiophosphorsäure-Derivat';
        187 :  f := 'Thiophosphorsäure';
        188 :  f := 'Thiophosphorsäureester';
        189 :  f := 'Thiophosphorsäurehalogenid';
        190 :  f := 'Thiophosphorsäureamid';
        191 :  f := 'Phosphonsäure-Derivat';
        192 :  f := 'Phosphonsäure';
        193 :  f := 'Phosphonsäureester';
        194 :  f := 'Phosphin';
        195 :  f := 'Phosphinoxid';
        196 :  f := 'Boronsäure-Derivat';
        197 :  f := 'Boronsäure';
        198 :  f := 'Boronsäureester';
        199 :  f := 'Alken';
        200 :  f := 'Alkin';
        201 :  f := 'aromatische Verbindung';
        202 :  f := 'heterocyclische Verbindung';
        203 :  f := 'alpha-Aminosäure';
        204 :  f := 'alpha-Hydroxysäure';
      end;  // case
    end;   // if lang = lang_de
  mkfglabel := f;
end;

procedure init_pt;  // v0.4, v0.4d
begin  // fills the periodic table with element symbols; index = atomic number (Z)
  pt[1].el := 'H '; pt[2].el := 'He';  pt[3].el := 'Li';  pt[4].el := 'Be';  pt[5].el := 'B '; 
  pt[1].am := 1.0079; pt[2].am := 4.0026;  pt[3].am := 6.941;  pt[4].am := 9.01218;  pt[5].am := 10.81; 

  pt[6].el := 'C '; pt[7].el := 'N ';  pt[8].el := 'O ';  pt[9].el := 'F ';  pt[10].el := 'Ne'; 
  pt[6].am := 12.011; pt[7].am := 14.00674;  pt[8].am := 15.9994;  pt[9].am := 18.999840;  pt[10].am := 20.179; 

  pt[11].el := 'Na'; pt[12].el := 'Mg';  pt[13].el := 'Al';  pt[14].el := 'Si';  pt[15].el := 'P '; 
  pt[11].am := 22.98977; pt[12].am := 24.305;  pt[13].am := 26.98154;  pt[14].am := 28.086;  pt[15].am := 30.97376; 

  pt[16].el := 'S '; pt[17].el := 'Cl';  pt[18].el := 'Ar';  pt[19].el := 'K ';  pt[20].el := 'Ca'; 
  pt[16].am := 32.06; pt[17].am := 35.453;  pt[18].am := 39.948;  pt[19].am := 39.098;  pt[20].am := 40.08; 

  pt[21].el := 'Sc'; pt[22].el := 'Ti';  pt[23].el := 'V ';  pt[24].el := 'Cr';  pt[25].el := 'Mn'; 
  pt[21].am := 44.9559; pt[22].am := 47.90;  pt[23].am := 50.9414;  pt[24].am := 51.996;  pt[25].am := 54.9380; 

  pt[26].el := 'Fe'; pt[27].el := 'Co';  pt[28].el := 'Ni';  pt[29].el := 'Cu';  pt[30].el := 'Zn'; 
  pt[26].am := 55.847; pt[27].am := 58.9332;  pt[28].am := 58.70;  pt[29].am := 63.546;  pt[30].am := 65.38; 

  pt[31].el := 'Ga'; pt[32].el := 'Ge';  pt[33].el := 'As';  pt[34].el := 'Se';  pt[35].el := 'Br'; 
  pt[31].am := 	69.72; pt[32].am := 72.59;  pt[33].am := 74.9216;  pt[34].am := 78.96;  pt[35].am := 79.904; 

  pt[36].el := 'Kr'; pt[37].el := 'Rb';  pt[38].el := 'Sr';  pt[39].el := 'Y ';  pt[40].el := 'Zr'; 
  pt[36].am := 83.80; pt[37].am := 85.4678;  pt[38].am := 87.62;  pt[39].am := 88.9059;  pt[40].am := 91.22; 

  pt[41].el := 'Nb'; pt[42].el := 'Mo';  pt[43].el := 'Tc';  pt[44].el := 'Ru';  pt[45].el := 'Rh'; 
  pt[41].am := 92.9064; pt[42].am := 95.94;  pt[43].am := 97;  pt[44].am := 101.07;  pt[45].am := 102.9055; 

  pt[46].el := 'Pd'; pt[47].el := 'Ag';  pt[48].el := 'Cd';  pt[49].el := 'In';  pt[50].el := 'Sn'; 
  pt[46].am := 106.4; pt[47].am := 106.4;  pt[48].am := 112.40;  pt[49].am := 114.82;  pt[50].am := 118.69; 

  pt[51].el := 'Sb'; pt[52].el := 'Te';  pt[53].el := 'I ';  pt[54].el := 'Xe';  pt[55].el := 'Cs'; 
  pt[51].am := 121.75; pt[52].am := 127.60;  pt[53].am := 126.9045;  pt[54].am := 131.30;  pt[55].am := 132.9054; 

  pt[56].el := 'Ba'; pt[57].el := 'La';  pt[58].el := 'Ce';  pt[59].el := 'Pr';  pt[60].el := 'Nd'; 
  pt[56].am := 137.34; pt[57].am := 138.9055;  pt[58].am := 140.12;  pt[59].am := 140.9077;  pt[60].am := 144.24; 

  pt[61].el := 'Pm'; pt[62].el := 'Sm';  pt[63].el := 'Eu';  pt[64].el := 'Gd';  pt[65].el := 'Tb'; 
  pt[61].am := 145; pt[62].am := 150.4;  pt[63].am := 151.96;  pt[64].am := 157.25;  pt[65].am := 158.9254; 

  pt[66].el := 'Dy'; pt[57].el := 'Ho';  pt[68].el := 'Er';  pt[69].el := 'Tm';  pt[70].el := 'Yb'; 
  pt[66].am := 162.50; pt[57].am := 164.9304;  pt[68].am := 167.26;  pt[69].am := 168.9342;  pt[70].am := 173.04; 

  pt[71].el := 'Lu'; pt[72].el := 'Hf';  pt[73].el := 'Ta';  pt[74].el := 'W ';  pt[75].el := 'Re'; 
  pt[71].am := 174.96; pt[72].am := 178.49;  pt[73].am := 180.9479;  pt[74].am := 183.5;  pt[75].am := 186.207; 

  pt[76].el := 'Os'; pt[77].el := 'Ir';  pt[78].el := 'Pt';  pt[79].el := 'Au';  pt[80].el := 'Hg'; 
  pt[76].am := 190.2; pt[77].am := 192.22;  pt[78].am := 195.09;  pt[79].am := 196.9665;  pt[80].am := 200.59; 

  pt[81].el := 'Tl'; pt[82].el := 'Pb';  pt[83].el := 'Bi';  pt[84].el := 'Po';  pt[85].el := 'At'; 
  pt[81].am := 204.37; pt[82].am := 207.2;  pt[83].am := 208.9804;  pt[84].am := 209;  pt[85].am := 210; 

  pt[86].el := 'Rn'; pt[87].el := 'Fr';  pt[88].el := 'Ra';  pt[89].el := 'Ac';  pt[90].el := 'Th'; 
  pt[86].am := 222; pt[87].am := 223;  pt[88].am := 226.0254;  pt[89].am := 227;  pt[90].am := 232.0381; 

  pt[91].el := 'Pa'; pt[92].el := 'U ';  pt[93].el := 'Np';  pt[94].el := 'Pu';  pt[95].el := 'Am'; 
  pt[91].am := 231.0359; pt[92].am := 238.029;  pt[93].am := 237.0482;  pt[94].am := 244;  pt[95].am := 243; 

  pt[96].el := 'Cm'; pt[97].el := 'Bk';  pt[98].el := 'Cf';  pt[99].el := 'Es';  pt[100].el := 'Fm'; 
  pt[96].am := 247; pt[97].am := 247;  pt[98].am := 251;  pt[99].am := 254;  pt[100].am := 257; 

  pt[101].el := 'Md'; pt[102].el := 'No';  pt[103].el := 'Lr';  pt[104].el := 'Rf';  pt[105].el := 'Db'; 
  pt[101].am := 258; pt[102].am := 259;  pt[103].am := 262;  pt[104].am := 261.11;  pt[105].am := 	268; 

  pt[106].el := 'Sg'; pt[107].el := 'Bh';  pt[108].el := 'Hs';  pt[109].el := 'Mt';  pt[110].el := 'Ds'; 
  pt[106].am := 271; pt[107].am := 270;  pt[108].am := 269;  pt[109].am := 278;  pt[110].am := 281; 

  pt[111].el := 'Rg'; pt[112].el := 'Cn';  pt[113].el := 'Ut';  pt[114].el := 'Uq';  pt[115].el := 'Up'; 
  pt[111].am := 281; pt[112].am := 285;  pt[113].am := 286 ;  pt[114].am := 289;  pt[115].am := 289; 

  pt[116].el := 'Uh'; pt[117].el := 'Us';  pt[118].el := 'Uo';  // from 112 there should be 3 letters: Uub ...
  pt[116].am := 293; pt[117].am := 294;  pt[118].am := 294;  // from 112 there should be 3 letters: Uub ...
end;


function atomicnumber(q:str2):integer;  // v0.4
var
  i, res : integer;
begin
  res := 0;
  if (q <> '') then
    begin
      if length(q) < 2 then q := q + ' ';
      q[2] := lowercase(q[2]);  // v0.4e
      for i := 1 to max_atomicnum do
        begin
          if (q = pt[i].el) then res := i;
        end;
    end;
  if (q = 'Q ') or (q = 'A ') or (q = 'X ') then res := 999;
  atomicnumber := res;
end;


procedure init_globals;
var
  i : integer;
begin
  opt_verbose     := false;
  opt_debug       := false;
  opt_exact       := false;
  opt_stdin       := false;
  opt_text        := false;
  opt_code        := false;
  opt_bin         := false;
  opt_bitstring   := false;
  opt_molout      := false;
  opt_molstat     := false;
  opt_molstat_X   := false;
  opt_molstat_v   := false;  // new in v0.4d
  opt_xmdlout     := false;
  opt_strict      := false;  // new in v0.2f
  opt_metalrings  := false;  // new in v0.3
  opt_geom        := false;  // new in v0.3d
  opt_chiral      := false;  // new in v0.3f
  opt_fp          := false;  // new in v0.3m
  opt_iso         := false;  // new in v0.3p
  opt_chg         := false;  // new in v0.3p
  opt_rad         := false;  // new in v0.3p
  //cm_mdlmolfile   := false;
  found_arominfo  := false;
  found_querymol  := false;
  ndl_querymol    := false;
  opt_rs          := rs_sar;  // v0.3i
  //ringsearch_mode := rs_sar;
  rfile_is_open   := false;  // new in v0.2g
  ez_search       := false;  // new in v0.3d
  rs_search       := false;  // new in v0.3f
  ez_flag         := false;  // new in v0.3f
  chir_flag       := false;  // new in v0.3f
  rs_strict       := false;  // new in v0.3j
  n_Ctot := 0; n_Otot := 0; n_Ntot := 0;             // new in v0.3g
  ndl_n_Ctot := 0; ndl_n_Otot := 0; ndl_n_Ntot := 0; // new in v0.3g
  for i := 1 to max_fg do fg[i] := false;
  try
    getmem(molbuf,sizeof(molbuftype));
  except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end;
  ether_generic   := false;       // v0.3j
  amine_generic   := false;       // v0.3j
  hydroxy_generic := false;       // v0.3j
  fpformat        := fpf_decimal; // v0.3m
  fpindex         := 0;           // v0.3m
  fp_exacthit     := false;       // v0.3m
  fp_exactblock   := false;       // v0.3m
  tmfcode         := 0;           // v0.3m
  tmfmismatch     := false;       // v0.3m
  auto_ssr        := false;       // v0.3n
  keep_DT         := true;        // v0.3p
  opt_hfp         := false;       // v0.4
  hfpformat       := fpf_decimal; // v0.4
  opt_matchnum    := false;       // v0.4a
  opt_matchnum1   := false;       // v0.4a
  opt_pos         := false;       // v0.5
  use_gmm         := false;        // v0.4a
  valid_gmm       := false;        // v0.4a
  fglang          := -1;           // v0.5
  fgloc           := nil;          // v0.5
end;


procedure init_molstat(var mstat:molstat_rec);
begin
  (*
  with mstat do
    begin
      n_QA := 0; n_QB := 0; n_chg := 0;
      n_C1 := 0; n_C2 := 0; n_C  := 0;
      n_CHB1p := 0; n_CHB2p := 0; n_CHB3p := 0; n_CHB4 := 0;
      n_O2 := 0; n_O3  := 0;
      n_N1 := 0; n_N2 := 0; n_N3 := 0;
      n_S := 0; n_SeTe := 0;
      n_F := 0; n_Cl := 0; n_Br := 0; n_I := 0;
      n_P := 0; n_B := 0;
      n_Met := 0; n_X := 0;
      n_b1 := 0; n_b2 := 0; n_b3 := 0; n_bar := 0;
      n_C1O := 0; n_C2O := 0; n_CN := 0; n_XY := 0;
      n_r3 := 0; n_r4 := 0; n_r5 := 0; n_r6 := 0; n_r7 := 0;
      n_r8 := 0; n_r9 := 0; n_r10 := 0; n_r11 := 0; n_r12 := 0; n_r13p := 0;
      n_rN := 0; n_rN1 := 0; n_rN2 := 0; n_rN3p := 0;
      n_rO := 0; n_rO1 := 0; n_rO2p := 0;
      n_rS := 0; n_rX := 0;
      n_rAr := 0; n_rBz := 0; n_br2p := 0;
    end;
  *)
  fillchar(mstat,sizeof(molstat_rec),0);  // v0.3k
end;


procedure debugoutput(dstr:string);
begin
  if opt_debug then writeln(dstr);
end;


procedure left_trim(var trimstr:string);
begin
  while (length(trimstr)>0) and ((trimstr[1]=' ') or (trimstr[1]=TAB)) do
    delete(trimstr,1,1);
end;


function left_int(var trimstr:string):integer;
var
  numstr : string;
  auxstr : string;
  auxint, code : integer;
begin
  numstr := '-+0123456789';
  auxstr := '';
  auxint := 0;
  while (length(trimstr)>0) and ((trimstr[1]=' ') or (trimstr[1]=TAB)) do
    delete(trimstr,1,1);
  while (length(trimstr)>0) and (pos(trimstr[1],numstr)>0) do
    begin
      auxstr := auxstr + trimstr[1];
      delete(trimstr,1,1);
    end;
  val(auxstr,auxint,code);
  left_int := auxint;
end;


procedure clear_atom_tags;
var
  i : integer;
begin
  if n_atoms > 0 then
    for i := 1 to n_atoms do atom^[i].tag := false;
end;


procedure set_atom_tags;
var
  i : integer;
begin
  if n_atoms > 0 then
    for i := 1 to n_atoms do atom^[i].tag := true;
end;


procedure clear_ndl_atom_tags;
var
  i : integer;
begin
  if ndl_n_atoms > 0 then
    for i := 1 to ndl_n_atoms do ndl_atom^[i].tag := false;
end;


procedure set_ndl_atom_tags;
var
  i : integer;
begin
  if ndl_n_atoms > 0 then
    for i := 1 to ndl_n_atoms do ndl_atom^[i].tag := true;
end;


function count_tagged_ndl_heavyatoms:integer;
var
  i, n : integer;
begin
  n := 0;
  if (ndl_n_atoms > 0) then
    begin
      for i := 1 to ndl_n_atoms do
        begin
          if (ndl_atom^[i].heavy and ndl_atom^[i].tag) then inc(n);
        end;
    end;
  count_tagged_ndl_heavyatoms := n;
end;


procedure clear_gmm;  // v0.4a
begin
  fillchar(gmm^,sizeof(global_matchmatrix),false);
end;

procedure clear_gmm_total;  // v0.4b
begin
  fillchar(gmm_total^,sizeof(global_matchmatrix),false);
end;


procedure write_gmm;  // v0.4a
var
  i, j : integer;
  c : char;
begin
  if (n_atoms < 1) or (ndl_n_atoms < 1) then exit;
  for j := 1 to n_atoms do
    begin
      for i := 1 to ndl_n_atoms do
        begin
          if gmm^[i,j] then c := 'X' else c := '.';
          write(c);
        end;
        writeln;
    end;
end;

//============================= geometry functions ==========================

function dist3d(p1,p2:p_3d):double;
var
  res : double;
begin
  res    := sqrt(sqr(p1.x-p2.x) + sqr(p1.y-p2.y) + sqr(p1.z-p2.z));
  dist3d := res;
end;

(*
function is_cis(p1,p2,p3,p4:p_3d):boolean;  // new in v0.3d
var                         // just a simple, distance-based estimation
  total_dist  : double;     // instead of calculating the dihedral angle
  direct_dist : double;
  res         : boolean;
begin
  res := false;
  total_dist  := dist3d(p1,p2) + dist3d(p2,p3) + dist3d(p3,p4);
  direct_dist := dist3d(p1,p4);
  if (direct_dist < 0.78 * total_dist) then res := true;  // cutoff value of 0.78 was
  is_cis := res;                                          // experimentally determined
end;
*)
// function is_cis was replaced by a new one in v0.3h


function subtract_3d(p1,p2:p_3d):p_3d;
var
  p : p_3d;
begin
  p.x := p1.x - p2.x;
  p.y := p1.y - p2.y;
  p.z := p1.z - p2.z;
  subtract_3d := p;
end;


function add_3d(p1,p2:p_3d):p_3d;
var
  p : p_3d;
begin
  p.x := p1.x + p2.x;
  p.y := p1.y + p2.y;
  p.z := p1.z + p2.z;
  add_3d := p;
end;


procedure vec2origin(var p1,p2:p_3d);
var
  p : p_3d;
begin
  p := subtract_3d(p2,p1);
  p2 := p;
  p1.x := 0; p1.y := 0; p1.z := 0;
end;


function scalar_prod(p1,p2,p3:p_3d):double;
var
  p : p_3d;
  res : double;
begin
  p := subtract_3d(p2,p1);
  p2 := p;
  p := subtract_3d(p3,p1);
  p3 := p;
  p1.x := 0; p1.y := 0; p1.z := 0;
  res := p2.x*p3.x + p2.y*p3.y + p2.z*p3.z;
  scalar_prod := res;
end;


function cross_prod(p1,p2,p3:p_3d):p_3d;
var
  p : p_3d;
  orig_p1 : p_3d;
begin
  orig_p1 := p1;
  p := subtract_3d(p2,p1);
  p2 := p;
  p := subtract_3d(p3,p1);
  p3 := p;
  p.x := p2.y*p3.z - p2.z*p3.y;
  p.y := p2.z*p3.x - p2.x*p3.z;
  p.z := p2.x*p3.y - p2.y*p3.x;
  cross_prod := add_3d(orig_p1,p);
end;


function angle_3d(p1,p2,p3:p_3d):double;
var
  lp1,lp2,lp3 : p_3d;
  p : p_3d;
  res : double;
  magn_1, magn_2 : double;
  cos_phi : double;
begin
  res := 0;
  lp1 := p1; lp2 := p2; lp3 := p3;
  p := subtract_3d(lp2,lp1);
  lp2 := p;
  p := subtract_3d(lp3,lp1);
  lp3 := p;
  lp1.x := 0; lp1.y := 0; lp1.z := 0;
  magn_1 := dist3d(lp1,lp2);
  magn_2 := dist3d(lp1,lp3);
  if (magn_1 * magn_2 = 0) then
    begin   // emergency exit
      angle_3d := pi;
      exit;
    end;
  cos_phi := scalar_prod(lp1,lp2,lp3) / (magn_1 * magn_2);
  if cos_phi < -1 then cos_phi := -1;
  if cos_phi > 1  then cos_phi := 1;
  res := arccos(cos_phi);
  angle_3d := res;
end;


function torsion(p1,p2,p3,p4:p_3d):double;
var
  lp1,lp2,lp3,lp4 : p_3d;
  d1 : p_3d;
  c1,c2 : p_3d;
  res : double;
  c1xc2, c2xc1 : p_3d;
  dist1,dist2 : double;
  sign : double;
begin
  // copy everything into local variables
  lp1 := p1; lp2 := p2; lp3 := p3; lp4 := p4;
  // get the vector between the two central atoms
  d1 := subtract_3d(p3,p2);
  // shift the first atom parallel to be attached to p3 instead of p2
  lp1 := add_3d(p1,d1);
  // now get the cross product vectors
  c1 := cross_prod(lp3,lp2,lp1);
  c2 := cross_prod(lp3,lp2,lp4);
  res := angle_3d(p3,c1,c2);
  //now check if it is clockwise or anticlockwise:
  //first, make the cross products of the two cross products c1 and c2 (both ways)
  c1xc2 := cross_prod(lp3,c1,c2);
  c2xc1 := cross_prod(lp3,c2,c1);
  //next, get the distances from these points to our refernce point lp2
  dist1 := dist3d(lp2,c1xc2);
  dist2 := dist3d(lp2,c2xc1);
  if (dist1 <= dist2) then sign := 1 else sign := -1;
  torsion := sign*res;
end;


function ctorsion(p1,p2,p3,p4:p_3d):double;
// calculates "pseudo-torsion" defined by atoms 3 and 4, being both
// attached to atom 2, with respect to axis of atoms 1 and 2
var
  lp1,lp2,lp3,lp4 : p_3d;
  //d1 : p_3d;
  c1,c2 : p_3d;
  res : double;
  c1xc2, c2xc1 : p_3d;
  dist1,dist2 : double;
  sign : double;
begin
  // copy everything into local variables
  lp1 := p1; lp2 := p2; lp3 := p3; lp4 := p4;
  // get the cross product vectors
  c1 := cross_prod(lp2,lp1,lp3);
  c2 := cross_prod(lp2,lp1,lp4);
  res := angle_3d(p2,c1,c2);
  //now check if it is clockwise or anticlockwise:
  //first, make the cross products of the two cross products c1 and c2 (both ways)
  c1xc2 := cross_prod(lp2,c1,c2);
  c2xc1 := cross_prod(lp2,c2,c1);
  //next, get the distances from these points to our refernce point lp1
  dist1 := dist3d(lp1,c1xc2);
  dist2 := dist3d(lp1,c2xc1);
  if (dist1 <= dist2) then sign := 1 else sign := -1;
  ctorsion := sign*res;
end;


function is_cis(p1,p2,p3,p4:p_3d):boolean;  // new in v0.3h, uses the dihedral angle
var
  phi : double;
  res : boolean;
begin
  res := false;
  phi := torsion(p1,p2,p3,p4);
  if (abs(phi) < pi/2) then res := true;
  is_cis := res;
end;


//====================== end of geometry functions ==========================

procedure show_usage;
begin
  if progmode = pmMatchMol then
    begin
      writeln('matchmol version ',version,'  N. Haider, University of Vienna, 2003-2018');
      writeln('Usage: matchmol [options] <needle> <haystack>');
      writeln(' where <needle> and <haystack> are the two molecules to compare');
      writeln(' (supported formats: MDL *.mol or *.sdf, Alchemy *.mol, Sybyl *.mol2)');
      writeln(' options can be:');
      writeln('    -v  verbose output');
      writeln('    -x  exact match');
      writeln('    -s  strict comparison of atom and bond types (including ring check)');  // new in v0.2f, v0.3d
      writeln('    -r  force SSR (set of small rings) ring search mode');
      writeln('    -m  write matching molecule as MDL molfile to standard output');
      writeln('        (default output: record number + ":T" for hit  or ":F" for miss');
      writeln('    -M  accept metal atoms as ring members');
      writeln('    -n  additional output of atom numbers for matching atom pairs');
      writeln('    -N  like -n, but only for the first matching substructure found');
      writeln('    -g  check geometry of double bonds (E/Z)');
      writeln('    -G  check geometry of chiral centers (R/S)');
      writeln('    -C  ignore all chirality checks even if "chiral" flag is set in molfile');
      writeln('    -a  check charges strictly');         // 0.3p
      writeln('    -i  check isotopes strictly');        // 0.3p
      writeln('    -d  check radicals strictly');        // 0.3p*/
      writeln('    -f  fingerprint mode (1 haystack, multiple needles) with boolean output');
      writeln('    -F  fingerprint mode (1 haystack, multiple needles) with decimal output');
    end else
    begin
      writeln('checkmol version ',version,'  N. Haider, University of Vienna, 2003-2018');
      writeln('Usage: checkmol [options] <filename>');
      writeln(' where options can be:');
      writeln('    -l  print a list of fingerprint codes + explanation and exit');
      writeln('    -v  verbose output');
      writeln('    -r  force SSR (set of small rings) ring search mode');
      writeln('    -M  accept metal atoms as ring members');
      writeln('  and one of the following:');
      writeln('    -e  english text (common name of functional group; default)');
      writeln('    -d  german text (common name of functional group)');
      writeln('    -c  code (acronym-like code for functional group)');
      writeln('    -b  bitstring (in decimal format) representing the presence of each group');
      writeln('    -s  the ASCII representation of the above bitstring, i.e. 0s and 1s)');
      writeln('    -p  lists the position of each functional group (atom number of key atom)');
      writeln('    -x  print molecular statistics (number of various atom types, bond types,');
      writeln('        ring sizes, etc.');
      writeln('    -X  same as above, listing all records (even if 0) as comma-separated list');
      writeln('    -a  count charges in fingerprint');  // 0.3p
      writeln('    -m  write MDL molfile (with special encoding for aromatic atoms/bonds)');
      writeln('    -h  hashed fingerprint mode with boolean output');
      writeln('    -H  hashed fingerprint mode with decimal output');
      writeln(' options can be combined like -vc');
      writeln(' <filename> specifies any file in the formats supported');
      writeln(' (MDL *.mol, Alchemy *.mol, Sybyl *.mol2), the filename "-" (without quotes)');
      writeln(' specifies standard input');
    end;
    // the "debug" option (-D) remains undocumented
end;


procedure list_molstat_codes;
begin
  writeln('n_atoms:     number of heavy atoms');
  writeln('n_bonds:     number of bonds between non-H atoms');
  writeln('n_rings:     number of rings');
  writeln('n_QA:        number of query atoms');
  writeln('n_QB:        number of query bonds');
  writeln('n_chg:       number of charges');
  writeln('n_C1:        number of sp-hybridized carbon atoms');
  writeln('n_C2:        number of sp2-hybridized carbon atoms');
  writeln('n_C:         total number of carbon atoms');
  writeln('n_CHB1p:     number of carbon atoms with at least 1 bond to a hetero atom');
  writeln('n_CHB2p:     number of carbon atoms with at least 2 bonds to a hetero atom');
  writeln('n_CHB3p:     number of carbon atoms with at least 3 bonds to a hetero atom');
  writeln('n_CHB4:      number of carbon atoms with 4 bonds to a hetero atom');
  writeln('n_O2:        number of sp2-hybridized oxygen atoms');
  writeln('n_O3:        number of sp3-hybridized oxygen atoms');
  writeln('n_N1:        number of sp-hybridized nitrogen atoms');
  writeln('n_N2:        number of sp2-hybridized nitrogen atoms');
  writeln('n_N3:        number of sp3-hybridized nitrogen atoms');
  writeln('n_S:         number of sulfur atoms');
  writeln('n_SeTe:      total number of selenium and tellurium atoms');
  writeln('n_F:         number of fluorine atoms');
  writeln('n_Cl:        number of chlorine atoms');
  writeln('n_Br:        number of bromine atoms');
  writeln('n_I:         number of iodine atoms');
  writeln('n_P:         number of phosphorus atoms');
  writeln('n_B:         number of boron atoms');
  writeln('n_Met:       total number of metal atoms');
  writeln('n_X:         total number of "other" atoms (not listed above) and halogens');
  writeln('n_b1:        number of single bonds');
  writeln('n_b2:        number of double bonds');
  writeln('n_b3:        number of triple bonds');
  writeln('n_bar:       number of aromatic bonds');
  writeln('n_C1O:       number of C-O single bonds');
  writeln('n_C2O:       number of C=O double bonds');
  writeln('n_CN:        number of C/N bonds (any type)');
  writeln('n_XY:        number of heteroatom/heteroatom bonds (any type)');
  writeln('n_r3:        number of 3-membered rings');
  writeln('n_r4:        number of 4-membered rings');
  writeln('n_r5:        number of 5-membered rings');
  writeln('n_r6:        number of 6-membered rings');
  writeln('n_r7:        number of 7-membered rings');
  writeln('n_r8:        number of 8-membered rings');
  writeln('n_r9:        number of 9-membered rings');
  writeln('n_r10:       number of 10-membered rings');
  writeln('n_r11:       number of 11-membered rings');
  writeln('n_r12:       number of 12-membered rings');
  writeln('n_r13p:      number of 13-membered or larger rings');
  writeln('n_rN:        number of rings containing nitrogen (any number)');
  writeln('n_rN1:       number of rings containing 1 nitrogen atom');
  writeln('n_rN2:       number of rings containing 2 nitrogen atoms');
  writeln('n_rN3p:      number of rings containing 3 or more nitrogen atoms');
  writeln('n_rO:        number of rings containing oxygen (any number)');
  writeln('n_rO1:       number of rings containing 1 oxygen atom');
  writeln('n_rO2p:      number of rings containing 2 or more oxygen atoms');
  writeln('n_rS:        number of rings containing sulfur (any number)');
  writeln('n_rX:        number of heterocycles (any type)');
  writeln('n_rar:       number of aromatic rings (any type)');
  {$IFDEF extended_molstat}
  writeln('n_rbz:       number of benzene rings');
  writeln('n_br2p:      number of bonds belonging to two or more rings');
  writeln('n_psg01:     number of atoms belonging to group 1 of the periodic system');
  writeln('n_psg02:     number of atoms belonging to group 2 of the periodic system');
  writeln('n_psg13:     number of atoms belonging to group 13 of the periodic system');
  writeln('n_psg14:     number of atoms belonging to group 14 of the periodic system');
  writeln('n_psg15:     number of atoms belonging to group 15 of the periodic system');
  writeln('n_psg16:     number of atoms belonging to group 16 of the periodic system');
  writeln('n_psg17:     number of atoms belonging to group 17 of the periodic system');
  writeln('n_psg18:     number of atoms belonging to group 18 of the periodic system');
  writeln('n_pstm:      number of atoms belonging to the transition metals');
  writeln('n_psla:      number of atoms belonging to the lanthanides or actinides');
  {$ENDIF}
end;


procedure parse_args;
var
  p : integer;
  parstr : string;
  tmpstr : string;
  l : integer;
begin
  tmpstr := '';
  opt_none := true;
  if progmode = pmCheckMol then
    begin
      for p := 1 to paramcount do
        begin
          parstr := paramstr(p);
          if (parstr = '-l') then   // new in v0.3l
            begin
              list_molstat_codes;
              halt(0);
            end;
          if p < paramcount then
            begin
              if (pos('-',parstr)=1) and (p < (paramcount)) then
                begin
                  tmpstr := paramstr(p);
                  left_trim(tmpstr);
                  l := 0;
                  if pos('v',tmpstr) > 0 then inc(l);
                  if pos('D',tmpstr) > 0 then inc(l);
                  if pos('r',tmpstr) > 0 then inc(l);
                  if pos('a',tmpstr) > 0 then inc(l);  // v0.3p
                  if pos('M',tmpstr) > 0 then inc(l);  // new in v0.3
                  if pos('H',tmpstr) > 0 then inc(l);  // new in v0.4
                  if pos('x',tmpstr) > 0 then inc(l);  // v0.4
                  if pos('X',tmpstr) > 0 then inc(l);  // v0.4
                  if pos('b',tmpstr) > 0 then inc(l);  // v0.4
                  if pos('s',tmpstr) > 0 then inc(l);  // v0.4
                  if pos('c',tmpstr) > 0 then inc(l);  // v0.4
                  if pos('p',tmpstr) > 0 then inc(l);  // v0.5
                  if (length(tmpstr) > (2 + l)) then
                     begin
                       show_usage;
                       halt(1);
                     end;
                  opt_none := false;
                  if pos('M',tmpstr)>0 then opt_metalrings  := true;
                  if pos('v',tmpstr)>0 then opt_verbose     := true;
                  {$IFDEF debug}
                  if pos('D',tmpstr)>0 then opt_debug       := true;
                  {$ENDIF}
                  if pos('e',tmpstr)>0 then 
                    begin
                      opt_text        := true;
                      fglang := lang_en;  // v0.5
                      if pos('p',tmpstr)>0 then 
                        begin
                          opt_pos := true;
                          opt_text := false;
                        end;
                    end else
                    begin   // not opt_text
                      if pos('d',tmpstr)>0 then 
                        begin
                          opt_text_de := true;
                          fglang := lang_de;
                          if pos('p',tmpstr)>0 then 
                            begin
                              opt_pos := true;
                              opt_text_de := false;
                            end;
                        end else
                        begin   // not opt_text_de
                          if pos('p',tmpstr)>0 then 
                            begin
                              opt_pos := true;
                              fglang := -1;
                            end else
                            begin   // not opt_text_de
                              if pos('c',tmpstr)>0 then opt_code    := true else
                                begin
                                  if pos('b',tmpstr)>0 then opt_bin := true else
                                  if pos('s',tmpstr)>0 then opt_bitstring := true;
                                end;
                              if pos('x',tmpstr)>0 then 
                                begin
                                  opt_molstat := true;
                                  if pos('xx',tmpstr)>0 then opt_molstat_v := true;  // v0.4d
                                end;
                              if pos('r',tmpstr)>0 then opt_rs      := rs_ssr;
                              if pos('a',tmpstr)>0 then opt_chg     := true;  // v0.3p
                              if pos('X',tmpstr)>0 then
                                begin
                                  opt_molstat   := true;
                                  opt_molstat_X := true;
                                end;
                              if pos('h',parstr)>0 then // new in v0.4
                                begin
                                  opt_hfp  := true;
                                  hfpformat := fpf_boolean;
                                end;
                              if pos('H',parstr)>0 then // new in v0.4
                                begin
                                  opt_hfp  := true;
                                  hfpformat := fpf_decimal;
                                end;
                              if pos('m',tmpstr)>0 then
                                begin
                                  opt_text      := false;
                                  opt_text_de   := false;
                                  opt_bin       := false;
                                  opt_bitstring := false;
                                  opt_code      := false;
                                  opt_pos       := false;  // v0.5
                                  opt_molstat   := false;
                                  opt_hfp       := false;
                                  opt_xmdlout   := true;
                                end;
                            end;  // not opt_pos
                        end;   // not opt_text_de
                    end;       // not opt_text
                  molfilename := tmpstr;
                end;
            end else
            begin
              if (pos('-',parstr)=1) then
                begin
                  if length(parstr)>1 then
                    begin
                      show_usage;
                      halt(1);
                    end else
                    begin
                      opt_stdin := true;
                    end;
                end else
                  begin
                    opt_stdin := false;
                    molfilename := parstr;
                  end;
            end;
        end;  // for p := 1 to paramcount
      if (opt_text = false) and (opt_text_de = false) and (opt_code = false) and
         (opt_bin = false) and  (opt_bitstring = false) and (opt_molstat = false) and
         (opt_molstat_X = false) and (opt_xmdlout = false) and (opt_chg = false) and
         (opt_hfp = false) and (opt_pos = false)    // v0.5
        then opt_none := true;
      if opt_pos then  // v0.5
        begin
          if opt_text    then fglang := lang_en;
          if opt_text_de then fglang := lang_de;
          opt_text       := false;
          opt_text_de    := false;
          opt_code       := false;
          opt_bin        := false;
          opt_bitstring  := false;
          opt_molstat    := false;
          opt_molstat_X  := false;
          opt_xmdlout    := false;
          opt_hfp        := false;
        end;
    end;
  if progmode = pmMatchMol then
    begin
      ndl_molfilename := '';
      molfilename := '';
      for p := 1 to paramcount do
        begin
          parstr := paramstr(p);
          if (p = 1) or (p < paramcount) then   // v0.4e
            begin
              if (pos('-',parstr)=1) then
                begin
                  if pos('v',parstr)>1 then opt_verbose     := true;
                  {$IFDEF debug}
                  if pos('D',parstr)>1 then opt_debug       := true;
                  {$ENDIF}
                  if pos('x',parstr)>1 then opt_exact       := true;
                  if pos('s',parstr)>1 then opt_strict      := true;  // new in v0.2f
                  if pos('m',parstr)>1 then opt_molout      := true;
                  if pos('n',parstr)>1 then 
                    begin
                      use_gmm         := true;  // new in v0.4a
                      opt_matchnum    := true;  // new in v0.4a
                      opt_matchnum1   := false;  // new in v0.4a
                    end;
                  if pos('N',parstr)>1 then 
                    begin
                      use_gmm         := true;  // new in v0.4a
                      opt_matchnum1   := true;  // new in v0.4a
                      opt_matchnum    := false;  // new in v0.4a
                    end;
                  if pos('r',parstr)>1 then opt_rs          := rs_ssr;
                  if pos('a',parstr)>1 then opt_chg         := true;  // v0.3p
                  if pos('i',parstr)>1 then opt_iso         := true;  // v0.3p
                  if pos('d',parstr)>1 then opt_rad         := true;  // v0.3p
                  if pos('M',parstr)>0 then opt_metalrings  := true;  // new in v0.3
                  if pos('g',parstr)>0 then opt_geom        := true;  // new in v0.3d
                  if pos('G',parstr)>0 then opt_chiral      := true;  // new in v0.3f
                  if pos('C',parstr)>0 then opt_nochirality := true;  // new in v0.4e
                  if pos('f',parstr)>0 then // new in v0.3m
                    begin
                      opt_fp   := true;
                      fpformat := fpf_boolean;
                    end;
                  if pos('F',parstr)>0 then // new in v0.3m
                    begin
                      opt_fp   := true;
                      fpformat := fpf_decimal;
                    end;
                  if pos('h',parstr)>1 then
                    begin
                      show_usage;
                      halt(0);
                    end;
                end else ndl_molfilename := parstr;
              end;
          if p = (paramcount - 1) then
            begin
              if (pos('-',parstr) <> 1) then ndl_molfilename := parstr;
            end;
          if p = paramcount then
            begin
              if parstr <> '-' then molfilename := parstr else opt_stdin := true;
            end;
        end;
      if opt_nochirality then opt_chiral := false;  // v0.4e
      if opt_geom then ez_search := true;  // v0.3d
      if opt_chiral then rs_search := true;  // v0.3f
      if (opt_chiral and opt_strict and (opt_exact or opt_fp)) then 
        rs_strict := true; // new in v0.3j, v0.3m
      if opt_strict then use_gmm := true;  // v0.4b
      if opt_fp then     // v0.3m
        begin
          opt_molout := false;
          opt_exact  := false;
          opt_matchnum  := false;  // v0.4a
          opt_matchnum1 := false;  // v0.4a
          use_gmm       := false;  // v0.4a
        end;
    end;  // progmode = pmMatchMol
  ringsearch_mode := opt_rs;  // v0.3i
end;


//============== input-related functions & procedures =====================

function get_filetype(f:string):string;
var
  rline : string;
  auxstr : string;
  i : integer;
  mdl1 : boolean;
  ri : integer;
  sepcount : integer;
begin
  auxstr := 'unknown';
  i := li; mdl1 := false;
  ri := li -1;
  sepcount := 0;
  while (ri < molbufindex) and (sepcount < 1) do
    begin
      inc(ri);
      rline := molbuf^[ri];
      if (pos('$$$$',rline)>0) then inc(sepcount);
      if (i = li) and (copy(rline,7,5)='ATOMS')
                 and (copy(rline,20,5)='BONDS')
                 and (copy(rline,33,7)='CHARGES') then
        begin
          auxstr := 'alchemy';
        end;
      if (i = li+3) // and (copy(rline,31,3)='999')
                 and (copy(rline,35,5)='V2000')      then mdl1 := true;
      if (i = li+1) and (copy(rline,3,6)='-ISIS-')      then mdl1 := true;
      if (i = li+1) and (copy(rline,3,8)='WLViewer')    then mdl1 := true;
      if (i = li+1) and (copy(rline,3,8)='CheckMol')    then mdl1 := true;
      if (i = li+1) and (copy(rline,3,8)='CATALYST') then
        begin
          mdl1 := true;
          auxstr := 'mdl';
        end;
      if (pos('M  END',rline)=1) or mdl1 then
        begin
          auxstr := 'mdl';
        end;
      if pos('@<TRIPOS>MOLECULE',rline)>0 then
        begin
          auxstr := 'sybyl';
        end;
      inc(i);
    end;
  // new in v0.2j: try to identify non-conformant SD-files
  if (auxstr = 'unknown') and (sepcount > 0) then auxstr := 'mdl';
  get_filetype := auxstr;
end;


procedure zap_molecule;
begin
  try
    if atom <> nil then
      begin
        freemem(atom,n_atoms*sizeof(atom_rec));
        atom := nil;  // added in v0.3j
      end;
    if bond <> nil then 
      begin
        freemem(bond,n_bonds*sizeof(bond_rec));
        bond := nil;  // added in v0.3j
      end;
    if ring <> nil then 
      begin
        freemem(ring,sizeof(ringlist));
        ring := nil;  // added in v0.3j
      end;
    if ringprop <> nil then 
      begin
        freemem(ringprop,sizeof(ringprop_type));
        ringprop := nil;  // added in v0.3j
      end;
    if fgloc <> nil then
      begin
        freemem(fgloc,sizeof(fgloctype));
        fgloc := nil;  // added in v0.5
      end;
  except
    on e:Einvalidpointer do begin end;
  end;
  n_atoms := 0;
  n_bonds := 0;
  n_rings := 0;
end;


procedure zap_needle;
begin
  try
    if ndl_atom <> nil then 
      begin
        freemem(ndl_atom,ndl_n_atoms*sizeof(atom_rec));
        ndl_atom := nil;  // added in v0.3j
      end;
    if ndl_bond <> nil then 
      begin
        freemem(ndl_bond,ndl_n_bonds*sizeof(bond_rec));
        ndl_bond := nil;  // added in v0.3j
      end;
    if ndl_ring <> nil then 
      begin
        freemem(ndl_ring,sizeof(ringlist));
        ndl_ring := nil;  // added in v0.3j
      end;
    if ndl_ringprop <> nil then 
      begin
        freemem(ndl_ringprop,sizeof(ringprop_type)); // fixed in v0.3g
        ndl_ringprop := nil;  // added in v0.3j
      end;
  except
    on e:Einvalidpointer do begin end;
  end;
  ndl_n_atoms := 0;
  ndl_n_bonds := 0;
  ndl_n_rings := 0;
end;


procedure zap_tmp;
begin
  try
    if tmp_atom <> nil then 
      begin
        freemem(tmp_atom,tmp_n_atoms*sizeof(atom_rec));
        tmp_atom := nil;  // added in v0.3j
      end;
    if tmp_bond <> nil then 
      begin
        freemem(tmp_bond,tmp_n_bonds*sizeof(bond_rec));
        tmp_bond := nil;  // added in v0.3j
      end;
    if tmp_ring <> nil then 
      begin
        freemem(tmp_ring,sizeof(ringlist));
        tmp_ring := nil;  // added in v0.3j
      end;
    if tmp_ringprop <> nil then 
      begin
        freemem(tmp_ringprop,sizeof(ringprop_type)); // fixed in v0.3g
        tmp_ringprop := nil;  // added in v0.3j
      end;
  except
    on e:Einvalidpointer do begin end;
  end;
  tmp_n_atoms := 0;
  tmp_n_bonds := 0;
  tmp_n_rings := 0;
end;


function is_heavyatom(id:integer):boolean;
var
  r  : boolean;
  el : str2;
begin
  r  := true;
  el := atom^[id].element;
  if (el = 'DU') or (el = 'LP') then r := false;
  (*
  if (progmode = pmCheckMol) and ((el = 'H ') and 
    (atom^[id].nucleon_number < 2)) then r := false; // v0.3p
  if (progmode = pmMatchMol) and (opt_iso = false) and ((el = 'H ') and
    (atom^[id].nucleon_number < 2)) then r := false; // v0.3p
  *)
  if (el = 'H ') then
    begin 
      if (opt_iso = false) then r := false else
        begin
          if (atom^[id].nucleon_number < 2) then r := false;
        end;
    end;
  is_heavyatom := r;  
  // note: deuterium is regarded as a heavy atom if -i option is used;
  // this applies only to matchmol, in checkmol D/T is always non-heavy
end;

function is_trueheavyatom(id:integer):boolean;   // v0.4b
var
  r  : boolean;
  el : str2;
begin
  r  := true;
  el := atom^[id].element;
  if (el = 'DU') or (el = 'LP') then r := false;
  if (el = 'H ') or (el = 'D ') or (el = 'T ') then r := false;
  is_trueheavyatom := r;  
end;

function ndl_alkene_C(ba:integer):boolean;   // new in v0.3f
var
  res : boolean;
  i, ba2 : integer;
begin
  res := false;
  if (ndl_n_atoms > 0) and (ndl_n_bonds > 0) then
    begin
      for i := 1 to ndl_n_bonds do
        begin
          if ((ndl_bond^[i].a1 = ba) or (ndl_bond^[i].a2 = ba)) then
            begin
              if (ndl_bond^[i].a1 = ba) then ba2 := ndl_bond^[i].a2 else
                  ba2 := ndl_bond^[i].a1;
              if (ndl_atom^[ba].atype ='C2 ') and (ndl_atom^[ba2].atype ='C2 ') and
                 (ndl_bond^[i].btype ='D') and (ndl_bond^[i].arom = false) then
                 res := true;
            end;
        end;
    end;
  ndl_alkene_C := res;  
end;


function is_metal(id:integer):boolean;
var
  r  : boolean;
  el : str2;
begin
  r  := false;
  el := atom^[id].element;
  if (el = 'LI') or (el = 'NA') or (el = 'K ') or (el = 'RB') or (el = 'CS') or
     (el = 'BE') or (el = 'MG') or (el = 'CA') or (el = 'SR') or (el = 'BA') or
     (el = 'TI') or (el = 'ZR') or (el = 'CR') or (el = 'MO') or (el = 'MN') or
     (el = 'FE') or (el = 'CO') or (el = 'NI') or (el = 'PD') or (el = 'PT') or
     (el = 'SN') or (el = 'CU') or (el = 'AG') or (el = 'AU') or (el = 'ZN') or 
     (el = 'CD') or (el = 'HG') or (el = 'AL') or (el = 'SN') or (el = 'PB') or 
     (el = 'SB') or (el = 'BI')                                   // etc. etc.
    then r := true;
  is_metal := r;
end;


function get_nvalences(a_el:str2):integer;  // changed name and position in v0.3m
// preliminary version; should be extended to element/atomtype
var
  res : integer;
begin
  res := 1;
  if a_el = 'H ' then res := 1;
  if a_el = 'D ' then res := 1;  // v0.3n
  if a_el = 'C ' then res := 4;
  if a_el = 'N ' then res := 3;
  if a_el = 'O ' then res := 2;
  if a_el = 'S ' then res := 2;
  if a_el = 'SE' then res := 2;
  if a_el = 'TE' then res := 2;
  if a_el = 'P ' then res := 3;
  if a_el = 'F ' then res := 1;
  if a_el = 'CL' then res := 1;
  if a_el = 'BR' then res := 1;
  if a_el = 'I ' then res := 1;
  if a_el = 'B ' then res := 3;
  if a_el = 'LI' then res := 1;
  if a_el = 'NA' then res := 1;
  if a_el = 'K ' then res := 1;
  if a_el = 'CA' then res := 2;
  if a_el = 'SR' then res := 2;
  if a_el = 'MG' then res := 2;
  if a_el = 'FE' then res := 3;
  if a_el = 'MN' then res := 2;
  if a_el = 'HG' then res := 2;
  if a_el = 'SI' then res := 4;
  if a_el = 'SN' then res := 4;
  if a_el = 'ZN' then res := 2;
  if a_el = 'CU' then res := 2;
  if a_el = 'A ' then res := 4;
  if a_el = 'Q ' then res := 4;
  get_nvalences := res;
end;


function convert_type(oldtype : str4):str3;
var
  i : integer;
  newtype : str3;
begin
  newtype := copy(oldtype,1,3);
  for i := 1 to 3 do newtype[i] := upcase(newtype[i]);
  if newtype[1] = '~' then newtype := 'VAL';
  If newtype[1] = '*' then newtype := 'STR';
  convert_type := newtype;
end;


function convert_sybtype(oldtype : str5):str3;
var
  newtype : str3;
begin
//  NewType := Copy(OldType,1,3);
//  For i := 1 To 3 Do NewType[i] := UpCase(NewType[i]);
//  If NewType[1] = '~' Then NewType := 'VAL';
//  If NewType[1] = '*' Then NewType := 'STR';
  newtype := 'DU ';
  if oldtype = 'H    ' then newtype := 'H  ';
  if oldtype = 'C.ar ' then newtype := 'CAR';
  if oldtype = 'C.2  ' then newtype := 'C2 ';
  if oldtype = 'C.3  ' then newtype := 'C3 ';
  if oldtype = 'C.1  ' then newtype := 'C1 ';
  if oldtype = 'O.2  ' then newtype := 'O2 ';
  if oldtype = 'O.3  ' then newtype := 'O3 ';
  if oldtype = 'O.co2' then newtype := 'O2 ';
  if oldtype = 'O.spc' then newtype := 'O3 ';
  if oldtype = 'O.t3p' then newtype := 'O3 ';
  if oldtype = 'N.1  ' then newtype := 'N1 ';
  if oldtype = 'N.2  ' then newtype := 'N2 ';
  if oldtype = 'N.3  ' then newtype := 'N3 ';
  if oldtype = 'N.pl3' then newtype := 'NPL';
  if oldtype = 'N.4  ' then newtype := 'N3+';
  if oldtype = 'N.am ' then newtype := 'NAM';
  if oldtype = 'N.ar ' then newtype := 'NAR';
  if oldtype = 'F    ' then newtype := 'F  ';
  if oldtype = 'Cl   ' then newtype := 'CL ';
  if oldtype = 'Br   ' then newtype := 'BR ';
  if oldtype = 'I    ' then newtype := 'I  ';
  if oldtype = 'Al   ' then newtype := 'AL ';
  if oldtype = 'ANY  ' then newtype := 'A  ';
  if oldtype = 'Ca   ' then newtype := 'CA ';
  if oldtype = 'Du   ' then newtype := 'DU ';
  if oldtype = 'Du.C ' then newtype := 'DU ';
  if oldtype = 'H.spc' then newtype := 'H  ';
  if oldtype = 'H.t3p' then newtype := 'H  ';
  if oldtype = 'HAL  ' then newtype := 'Cl ';
  if oldtype = 'HET  ' then newtype := 'Q  ';
  if oldtype = 'HEV  ' then newtype := 'DU ';
  if oldtype = 'K    ' then newtype := 'K  ';
  if oldtype = 'Li   ' then newtype := 'LI ';
  if oldtype = 'LP   ' then newtype := 'LP ';
  if oldtype = 'Na   ' then newtype := 'NA ';
  if oldtype = 'P.3  ' then newtype := 'P3 ';
  if oldtype = 'S.2  ' then newtype := 'S2 ';
  if oldtype = 'S.3  ' then newtype := 'S3 ';
  if oldtype = 'S.o  ' then newtype := 'SO ';
  if oldtype = 'S.o2 ' then newtype := 'SO2';
  if oldtype = 'Si   ' then newtype := 'SI ';
  if oldtype = 'P.4  ' then newtype := 'P4 ';
  convert_sybtype := newtype;
end;


function convert_MDLtype(oldtype : str3):str3;
var
  newtype : str3;
begin
//  NewType := Copy(OldType,1,3);
//  For i := 1 To 3 Do NewType[i] := UpCase(NewType[i]);
//  If NewType[1] = '~' Then NewType := 'VAL';
//  If NewType[1] = '*' Then NewType := 'STR';
  newtype := 'DU ';
  if oldtype = 'H  ' then newtype := 'H  ';
  if oldtype = 'C  ' then newtype := 'C3 ';
  if oldtype = 'O  ' then newtype := 'O2 ';
  if oldtype = 'N  ' then newtype := 'N3 ';
  if oldtype = 'F  ' then newtype := 'F  ';
  if oldtype = 'Cl ' then newtype := 'CL ';
  if oldtype = 'Br ' then newtype := 'BR ';
  if oldtype = 'I  ' then newtype := 'I  ';
  if oldtype = 'Al ' then newtype := 'AL ';
  if oldtype = 'ANY' then newtype := 'A  ';
  if oldtype = 'Ca ' then newtype := 'CA ';
  if oldtype = 'Du ' then newtype := 'DU ';
  if oldtype = 'K  ' then newtype := 'K  ';
  if oldtype = 'Li ' then newtype := 'LI ';
  if oldtype = 'LP ' then newtype := 'LP ';
  if oldtype = 'Na ' then newtype := 'NA ';
  if oldtype = 'P  ' then newtype := 'P3 ';
  if oldtype = 'S  ' then newtype := 'S3 ';
  if oldtype = 'Si ' then newtype := 'SI ';
  if oldtype = 'P  ' then newtype := 'P4 ';
  if oldtype = 'A  ' then newtype := 'A  ';
  if oldtype = 'Q  ' then newtype := 'Q  ';
  convert_MDLtype := NewType;
end;


function get_element(oldtype:str4):str2;
var
  elemstr : string;
begin
  if oldtype = 'H   ' then elemstr := 'H ';
  if oldtype = 'D   ' then elemstr := 'D ';  // v0.3n
  if oldtype = 'CAR ' then elemstr := 'C ';
  if oldtype = 'C2  ' then elemstr := 'C ';
  if oldtype = 'C3  ' then elemstr := 'C ';
  if oldtype = 'C1  ' then elemstr := 'C ';
  if oldtype = 'O2  ' then elemstr := 'O ';
  if oldtype = 'O3  ' then elemstr := 'O ';
  if oldtype = 'O2  ' then elemstr := 'O ';
  if oldtype = 'O3  ' then elemstr := 'O ';
  if oldtype = 'O3  ' then elemstr := 'O ';
  if oldtype = 'N1  ' then elemstr := 'N ';
  if oldtype = 'N2  ' then elemstr := 'N ';
  if oldtype = 'N3  ' then elemstr := 'N ';
  if oldtype = 'NPL ' then elemstr := 'N ';
  if oldtype = 'N3+ ' then elemstr := 'N ';
  if oldtype = 'NAM ' then elemstr := 'N ';
  if oldtype = 'NAR ' then elemstr := 'N ';
  if oldtype = 'F   ' then elemstr := 'F ';
  if oldtype = 'CL  ' then elemstr := 'CL';
  if oldtype = 'BR  ' then elemstr := 'BR';
  if oldtype = 'I   ' then elemstr := 'I ';
  if oldtype = 'AL  ' then elemstr := 'AL';
  if oldtype = 'DU  ' then elemstr := 'DU';
  if oldtype = 'CA  ' then elemstr := 'CA';
  if oldtype = 'DU  ' then elemstr := 'DU';
  if oldtype = 'Cl  ' then elemstr := 'CL';
  if oldtype = 'K   ' then elemstr := 'K ';
  if oldtype = 'LI  ' then elemstr := 'LI';
  if oldtype = 'LP  ' then elemstr := 'LP';
  if oldtype = 'NA  ' then elemstr := 'NA';
  if oldtype = 'P3  ' then elemstr := 'P ';
  if oldtype = 'S2  ' then elemstr := 'S ';
  if oldtype = 'S3  ' then elemstr := 'S ';
  if oldtype = 'SO  ' then elemstr := 'S ';
  if oldtype = 'SO2 ' then elemstr := 'S ';
  if oldtype = 'SI  ' then elemstr := 'SI';
  if oldtype = 'P4  ' then elemstr := 'P ';
  if oldtype = 'A   ' then elemstr := 'A ';
  if oldtype = 'Q   ' then elemstr := 'Q ';
  get_element := elemstr;
end;


function get_sybelement(oldtype:str5):str2;
var
  i : integer;
  elemstr : str2;
begin
  if pos('.',oldtype)<2 then
    begin
      elemstr := copy(oldtype,1,2);
    end else
    begin
      elemstr := copy(oldtype,1,pos('.',oldtype)-1);
      if length(elemstr)<2 then elemstr := elemstr+' ';
    end;
  for i := 1 to 2 do elemstr[i] := upcase(elemstr[i]);
  get_sybelement := elemstr;
end;


function get_MDLelement(oldtype:str3):str2;
var
  i : integer;
  elemstr : str2;
begin
  elemstr := copy(oldtype,1,2);
  for i := 1 to 2 do elemstr[i] := upcase(elemstr[i]);
  if elemstr[1] = '~' then elemstr := '??';
  If elemstr[1] = '*' then elemstr := '??';
  get_MDLelement := elemstr;
end;

function get_sep_label(instring:string):string;    // v0.4b, v0.4c
var
  x, y : integer;
  res : string;
begin
  res := '';
  x := pos('#',instring);
  if (x > 0) then
    begin
      y := length(instring);
      if (y > 70) then y := 70;
      res := copy(instring,x,(y-x)+1);
    end;
  get_sep_label := res;
end;


procedure read_molfile(mfilename:string);  // reads ALCHEMY mol files
var
  n, code : integer;
  rline, tmpstr : string;
  xstr, ystr, zstr, chgstr : string;
  xval, yval, zval, chgval : single;
  a1str, a2str, elemstr : string;
  a1val, a2val : integer;
  ri : integer;
begin
  if n_atoms > 0 then zap_molecule;
  ri := li;
  rline := molbuf^[ri];
  tmpstr := copy(rline,1,5);
  val(tmpstr,n_atoms,code);
  tmpstr := copy(rline,14,5);
  val(tmpstr,n_bonds,code);
  molname := copy(rline,42,length(rline)-42);
  try
    getmem(atom,n_atoms*sizeof(atom_rec));
    fillchar(atom^,n_atoms*sizeof(atom_rec),0); // blank out all atom records with 0 at once; v0.3l
    getmem(bond,n_bonds*sizeof(bond_rec));
    getmem(ring,sizeof(ringlist));
    getmem(ringprop,sizeof(ringprop_type));
  except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end;
  n_heavyatoms := 0;
  n_trueheavyatoms := 0;  // v0.4b
  n_heavybonds := 0;
  n_Ctot       := 0;  // v0.3g
  n_Otot       := 0;  // v0.3g
  n_Ntot       := 0;  // v0.3g
  for n := 1 to n_atoms do
    begin
      (*
      with atom^[n] do
        begin
          x := 0; y := 0; z := 0;  // v0.3g
          formal_charge  := 0;
          real_charge    := 0;
          Hexp           := 0;
          Htot           := 0;
          neighbor_count := 0;
          ring_count     := 0;
          arom           := false;
          q_arom         := false;
          stereo_care    := false;
          heavy          := false;
          metal          := false;
          tag            := false;
          nucleon_number := 0;
          radical_type   := 0;
        end;
      *)
      inc(ri);
      rline := molbuf^[ri];
      atomtype := copy(rline,7,4);
      atomtype := upcase(atomtype);  // fixed in v0.3f
      elemstr  := get_element(atomtype);
      if (elemstr = 'C ') then inc(n_Ctot);
      if (elemstr = 'O ') then inc(n_Otot);
      if (elemstr = 'N ') then inc(n_Ntot);
      newatomtype := convert_type(atomtype);
      xstr := copy(rline,14,7);
      ystr := copy(rline,23,7);
      zstr := copy(rline,32,7);
      chgstr := copy(rline,43,7);
      val(xstr,xval,code);
      val(ystr,yval,code);
      val(zstr,zval,code);
      val(chgstr,chgval,code);
      with atom^[n] do
        begin
          element := elemstr;
          atype := newatomtype;
          if (elemstr = 'D ') then   // v0.3p
            begin
              element := 'H ';
              nucleon_number := 2;
              atype := 'DU ';
            end;
          if (elemstr = 'T ') then   // v0.3p
            begin
              element := 'H ';
              nucleon_number := 3;
              atype := 'DU ';
            end;
          x := xval; y := yval; z := zval; real_charge := chgval;
          if (is_heavyatom(n)) then 
	    begin
	      inc(n_heavyatoms);
	      heavy := true;
	      if is_metal(n) then metal := true;
	      if (is_trueheavyatom(n)) then inc(n_trueheavyatoms);  // v0.4b
	    end;
	  nvalences := get_nvalences(element);  // v0.3m  
        end;  // with
    end;  // for
  for n := 1 to n_bonds do
    begin
      inc(ri);
      rline := molbuf^[ri];
      a1str := copy(rline,9,3);
      a2str := copy(rline,15,3);
      val(a1str,a1val,code);
      if code <> 0 then beep;
      val(a2str,a2val,code);
      if code <> 0 then beep;
      with bond^[n] do
        begin
          a1 := a1val; a2 := a2val; btype := rline[20];
          ring_count := 0; arom := false; q_arom := false; // v0.3p
          topo := btopo_any; stereo := bstereo_any;
          mdl_stereo := 0;  // v0.3n
        end;
      if atom^[a1val].heavy and atom^[a2val].heavy then inc(n_heavybonds);
    end;
  fillchar(ring^,sizeof(ringlist),0);
  for n := 1 to max_rings do  // new in v0.3
    begin
      ringprop^[n].size     := 0;
      ringprop^[n].arom     := false;
      ringprop^[n].envelope := false;
    end;
  li := ri + 1;
  keep_DT := true;  // v0.3p
end;


procedure read_mol2file(mfilename:string);  // reads SYBYL mol2 files
var
  n, code : integer;
  sybatomtype : string[5];
  tmpstr, rline : string;
  xstr, ystr, zstr, chgstr : string;
  xval, yval, zval, chgval : single;
  a1str, a2str, elemstr : string;
  a1val, a2val : integer;
  ri : integer;
begin
  if n_atoms > 0 then zap_molecule;
  rline := '';
  ri := li -1;
  while (ri < molbufindex) and (pos('@<TRIPOS>MOLECULE',rline)=0) do
    begin
      inc(ri);
      rline := molbuf^[ri];
    end;
  if ri < molbufindex then
    begin
      inc(ri);
      molname := molbuf^[ri];
    end;
  if ri < molbufindex then
    begin
      inc(ri);
      rline := molbuf^[ri];
    end;
  tmpstr := copy(rline,1,5);
  val(tmpstr,n_atoms,code);
  tmpstr := copy(rline,7,5);
  val(tmpstr,n_bonds,code);
  try
    getmem(atom,n_atoms*sizeof(atom_rec));
    fillchar(atom^,n_atoms*sizeof(atom_rec),0); // blank out all atom records with 0 at once; v0.3l
    getmem(bond,n_bonds*sizeof(bond_rec));
    getmem(ring,sizeof(ringlist));
    getmem(ringprop,sizeof(ringprop_type));
  except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end;
  n_heavyatoms := 0;
  n_trueheavyatoms := 0;  // v0.4b
  n_heavybonds := 0;
  n_Ctot       := 0;  // v0.3g
  n_Otot       := 0;  // v0.3g
  n_Ntot       := 0;  // v0.3g
  while ((ri < molbufindex) and (pos('@<TRIPOS>ATOM',rline)=0)) do
    begin
      inc(ri);
      rline := molbuf^[ri];
    end;
  for n := 1 to n_atoms do
  begin
    (*
    with atom^[n] do
      begin
        x := 0; y := 0; z := 0;  // v0.3g
        formal_charge  := 0;
        real_charge    := 0;
        Hexp           := 0;
        Htot           := 0;
        neighbor_count := 0;
        ring_count     := 0;
        arom           := false;
        q_arom         := false;
        stereo_care    := false;
        heavy          := false;
        metal          := false;
        tag            := false;
        nucleon_number := 0;
        radical_type   := 0;
      end;
    *)
    if (ri < molbufindex) then
      begin
        inc(ri);
        rline := molbuf^[ri];
      end;
    sybatomtype := copy(rline,48,5);
    elemstr     := get_sybelement(sybatomtype);
    if (elemstr = 'C ') then inc(n_Ctot);
    if (elemstr = 'O ') then inc(n_Otot);
    if (elemstr = 'N ') then inc(n_Ntot);
    newatomtype := convert_sybtype(sybatomtype);
    xstr := copy(rline,18,9);
    ystr := copy(rline,28,9);
    zstr := copy(rline,38,9);
    chgstr := copy(rline,70,9);
    val(xstr,xval,code);
    val(ystr,yval,code);
    val(zstr,zval,code);
    val(chgstr,chgval,code);
    with atom^[n] do
      begin
        element := elemstr;
        atype := newatomtype;
        if (elemstr = 'D ') then   // v0.3p
          begin
            element := 'H ';
            nucleon_number := 2;
            atype := 'DU ';
          end;
        if (elemstr = 'T ') then   // v0.3p
          begin
            element := 'H ';
            nucleon_number := 3;
            atype := 'DU ';
          end;
        x := xval; y := yval; z := zval; 
        formal_charge := round(chgval); real_charge := chgval; // v0.5a
        if (is_heavyatom(n)) then 
	  begin
	    inc(n_heavyatoms);
	    heavy := true;
	    if is_metal(n) then metal := true;
	    if (is_trueheavyatom(n)) then inc(n_trueheavyatoms);  // v0.4b
	  end;
        nvalences := get_nvalences(element);  // v0.3m  
      end;
  end;
  while ((ri < molbufindex) and (pos('@<TRIPOS>BOND',rline)=0)) do
    begin
      inc(ri);
      rline := molbuf^[ri];
    end;
  for n := 1 to n_bonds do
  begin
    if (ri < molbufindex) then
      begin
        inc(ri);
        rline := molbuf^[ri];
      end;
    a1str := copy(rline,9,3);
    a2str := copy(rline,14,3);
    val(a1str,a1val,code);
    if code <> 0 then writeln(rline, #7);
    val(a2str,a2val,code);
    if code <> 0 then writeln(rline,#7);
    with bond^[n] do
      begin
        a1 := a1val; a2 := a2val;
        if rline[18] = '1' then btype := 'S';
        if rline[18] = '2' then btype := 'D';
        if rline[18] = '3' then btype := 'T';
        if rline[18] = 'a' then 
          begin    // v0.5a  treat 'am' (amide) bonds as single
            if rline[19] = 'm' then btype := 'S' else btype := 'A';
          end;
        ring_count := 0; arom := false; q_arom := false;  // v0.3p
        topo := btopo_any; stereo := bstereo_any;
        mdl_stereo := 0;  // v0.3n
      end;
    if atom^[a1val].heavy and atom^[a2val].heavy then inc(n_heavybonds);
  end;
  fillchar(ring^,sizeof(ringlist),0);
  for n := 1 to max_rings do  // new in v0.3
    begin
      ringprop^[n].size     := 0;
      ringprop^[n].arom     := false;
      ringprop^[n].envelope := false;
    end;
  li := ri + 1;
  keep_DT := true;  // v0.3p
end;


procedure read_charges(chgstring:string);
var
  a_id, a_chg : integer;
  n_chrg : integer;
  // typical example: a molecule with 2 cations + 1 anion
  // M  CHG  3   8   1  10   1  11  -1
begin
  if (pos('M  CHG',chgstring)>0) then
    begin
      delete(chgstring,1,6);
      left_trim(chgstring);
      n_chrg := left_int(chgstring);  // this assignment must be kept also in non-debug mode!
      {$IFDEF debug}
      if (n_chrg = 0) then debugoutput('strange... M  CHG present, but no charges found');
      {$ENDIF}
      while (length(chgstring) > 0) do
        begin
          a_id  := left_int(chgstring);
          a_chg := left_int(chgstring);
          if (a_id <> 0) and (a_chg <> 0) then atom^[a_id].formal_charge := a_chg;
        end;
    end;
end;


procedure read_isotopes(isotopestring:string);
var
  a_id, a_nucleon_number : integer;
  n_isotopes : integer;
  // typical example: a molecule with 3 isotopes
  // M  ISO  3   8   15  10   13  11  17
begin
  if (pos('M  ISO',isotopestring) > 0) then
    begin
      delete(isotopestring,1,6);
      left_trim(isotopestring);
      n_isotopes := left_int(isotopestring);  // this assignment must be kept also in non-debug mode!
      {$IFDEF debug}
      if (n_isotopes = 0) then debugoutput('strange... M  ISO with nucleon_numer = 0');
      {$ENDIF}
      while (length(isotopestring) > 0) do
        begin
          a_id  := left_int(isotopestring);
          a_nucleon_number := left_int(isotopestring);
          if (a_id <> 0) and (a_nucleon_number > 0) then 
            begin
              atom^[a_id].nucleon_number := a_nucleon_number;
              if (atom^[a_id].element = 'H ') and (a_nucleon_number > 1) then
                begin
                  keep_DT := false;
                  if opt_iso then
                    begin
                      atom^[a_id].heavy := true;
                      inc(n_heavyatoms);
                    end;
                  atom^[a_id].atype := 'DU ';
                end;
            end;
        end;
    end;
end;


procedure read_radicals(radstring:string);
var
  a_id, a_rad : integer;
  n_rads : integer;
  // typical example: a molecule with a radical
  // M  RAD  1   8   2
begin
  if (pos('M  RAD',radstring) > 0) then
    begin
      delete(radstring,1,6);
      left_trim(radstring);
      n_rads := left_int(radstring);  // this assignment must be kept also in non-debug mode!
      {$IFDEF debug}
      if (n_rads = 0) then debugoutput('strange... M  RAD present, but no radicals found');
      {$ENDIF}
      while (length(radstring) > 0) do
        begin
          a_id  := left_int(radstring);
          a_rad := left_int(radstring);
          if (a_id <> 0) and (a_rad <> 0) then atom^[a_id].radical_type := a_rad;
        end;
    end;
end;


procedure read_MDLmolfile(mfilename:string);  // reads MDL mol files
var
  n, v, tmp_n_atoms, tmp_n_bonds, code : integer;  // v0.3l
  rline, tmpstr : string;
  xstr, ystr, zstr, chgstr : string;
  xval, yval, zval, chgval : single;
  a1str, a2str, elemstr : string;
  a1val, a2val : integer;
  ri,rc,bt,bs : integer;
  sepcount : integer;
  i : integer;              // v0.3j
  clearcharges : boolean;   // v0.3j
begin
  found_querymol := false;  // v0.3p
  clearcharges := true;     // v0.3j
  if n_atoms > 0 then zap_molecule;
  //cm_mdlmolfile := false;
  rline := '';
  ri := li;
  molname := molbuf^[ri];            // line 1
  if ri < molbufindex then inc(ri);  // line 2
  rline   := molbuf^[ri];
  if pos('CheckMol',rline)=3 then
    begin
      //cm_mdlmolfile := true;
      found_arominfo := true;
      tmfcode := 1;  // v0.3m (begin)
      code := 0;
      if (length(rline) >= 39) and (pos('TMF',rline) = 35) then  // v0.3m; encoding of tweaklevel
        begin
          tmpstr := copy(rline,38,2);
          val(tmpstr,tmfcode,code);
        end;
      if (code <> 0) or (tmfcode <> tweaklevel) then 
          tmfmismatch := true else tmfmismatch := false;
      if ((pos(':r0',rline) >= 40) and (ringsearch_mode <> rs_sar)) or
         ((pos(':r1',rline) >= 40) and (ringsearch_mode <> rs_ssr)) then tmfmismatch := true;
      if ((pos(':m0',rline) >= 40) and (opt_metalrings = true)) or
         ((pos(':m1',rline) >= 40) and (opt_metalrings = false)) then tmfmismatch := true;
      {$IFDEF debug}
      if tmfmismatch then 
        debugoutput('"tweaked" molfile: version mismatch') else
        debugoutput('"tweaked" molfile: version OK');
      {$ENDIF}
      // v0.3m (end)
    end;
  if ri < molbufindex then inc(ri);  // line 3
  rline   := molbuf^[ri];
  molcomment := rline;
  if ri < molbufindex then inc(ri);  // line 4
  rline := molbuf^[ri];
  tmpstr := copy(rline,1,3);
  val(tmpstr,n_atoms,code);
  tmpstr := copy(rline,4,3);
  val(tmpstr,n_bonds,code);
  tmpstr := copy(rline,10,3);   // if it is a CheckMol-tweaked molfile, this is the number of rings
  n_cmrings := 0;
  val(tmpstr,n_cmrings,code);
  if code <> 0 then n_cmrings := 0;
  // do some range checking for n_atoms, n_bonds; new in v0.3l
  tmp_n_atoms := n_atoms;
  if (n_atoms > max_atoms) then n_atoms := max_atoms;
  if (n_atoms < 0) then n_atoms := 0;
  tmp_n_bonds := n_bonds;
  if (n_bonds > max_bonds) then n_bonds := max_bonds;
  if (n_bonds < 0) then n_bonds := 0;
  if (n_atoms = 0) and opt_verbose then  // v0.3l
    begin
      writeln('WARNING: Possible NoStruct read!');
      writeln('NoStructs are proprietary, obsolete and dangerous.');
    end;
  try
    getmem(atom,n_atoms*sizeof(atom_rec));
    fillchar(atom^,n_atoms*sizeof(atom_rec),0); // blank out all atom records with 0 at once;
    getmem(bond,n_bonds*sizeof(bond_rec));      // this would be only one calloc() in C;  v0.3l
    getmem(ring,sizeof(ringlist));
    getmem(ringprop,sizeof(ringprop_type));
  except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        //close(molfile);
        halt(4);
        exit;
      end;
  end;
  // check for the chirality flag
  if ((length(rline) > 14) and (rline[15]='1')) then chir_flag := true;  // new in v0.3f
  n_heavyatoms := 0;
  n_trueheavyatoms := 0;  // v0.4b
  n_heavybonds := 0;
  n_Ctot       := 0;  // v0.3g
  n_Otot       := 0;  // v0.3g
  n_Ntot       := 0;  // v0.3g
  if (n_atoms > 0) then  // v0.3l
    begin
      for n := 1 to tmp_n_atoms do
        begin
          if (n <= max_atoms) then v := n else v := max_atoms;  // just for safety; v0.3l
          (*
          with atom^[v] do
            begin
              x := 0; y := 0; z := 0;  // v0.3g
              formal_charge  := 0;
              real_charge    := 0;
              Hexp           := 0;
              Htot           := 0;
              neighbor_count := 0;
              ring_count     := 0;
              arom           := false;
              q_arom         := false;
              stereo_care    := false;
              metal          := false;
              heavy          := false;
              tag            := false;
              nucleon_number := 0;
              radical_type   := 0;
            end;
          *) // replaced by fillchar() after getmem() (see above); v0.3l
          inc(ri);
          rline := molbuf^[ri];
          atomtype := copy(rline,32,3);
          elemstr  := get_MDLelement(atomtype);
          if (elemstr = 'C ') then inc(n_Ctot);
          if (elemstr = 'O ') then inc(n_Otot);
          if (elemstr = 'N ') then inc(n_Ntot);
          newatomtype := convert_MDLtype(atomtype);
          xstr := copy(rline,1,10);   // fixed in v0.3k (was: 2,9 etc.)
          ystr := copy(rline,11,10);
          zstr := copy(rline,21,10);
          //chgstr := '0';
          chgstr := copy(rline,37,3);   // new in v0.3j
          val(chgstr,chgval,code);
          if (chgval <> 0) then
            begin
              if (chgval >= 1) and (chgval <= 7) then
                chgval := 4 - chgval else chgval := 0;
            end;                        // end (v0.3j)
          val(xstr,xval,code);
          val(ystr,yval,code);
          val(zstr,zval,code);  // v0.3k: removed superfluous val(chgstr,chgval,code)
          with atom^[v] do
            begin
              element := elemstr;
              if (elemstr = 'A ') or (elemstr = 'Q ') or (elemstr = 'X ') 
                then found_querymol := true;  // 'X ' added in v0.3n
              atype := newatomtype;
              if (elemstr = 'D ') then   // v0.3p
                begin
                  keep_DT := true;
                  element := 'H ';
                  nucleon_number := 2;
                  atype := 'DU ';
                end;
              if (elemstr = 'T ') then   // v0.3p
                begin
                  keep_DT := true;
                  element := 'H ';
                  nucleon_number := 3;
                  atype := 'DU ';
                end;
              x := xval; y := yval; z := zval; 
              formal_charge := round(chgval); real_charge := 0; // v0.3j
              // read aromaticity flag from CheckMol-tweaked MDL molfile
              if (length(rline) > 37) and (rline[38] = '0') then
                begin
                  arom := true;
                  found_arominfo := true;
                end;
              // new in v0.3d: read stereo care flag
              if (length(rline) > 47) and (rline[48] = '1') then stereo_care := true;
              if (is_heavyatom(n)) then 
                begin
                  inc(n_heavyatoms);
                  heavy := true;
                  if is_metal(n) then metal := true;
                  if (is_trueheavyatom(n)) then inc(n_trueheavyatoms);  // v0.4b
                end;
              nvalences := get_nvalences(element);  // v0.3m                
            end;
        end;
    end;  // if (n_atoms > 0)...
  if (n_bonds > 0) then  // v0.3l
    begin
      for n := 1 to tmp_n_bonds do
        begin
          if (n <= max_bonds) then v := n else v := max_bonds;  // just for safety; v0.3l
          inc(ri);
          rline := molbuf^[ri];
          a1str := copy(rline,1,3);
          a2str := copy(rline,4,3);
          val(a1str,a1val,code);
          if code <> 0 then a1val := 1;  // v0.3l
          val(a2str,a2val,code);
          if code <> 0 then a2val := 1;  // v0.3l
          with bond^[v] do
            begin
              a1 := a1val; a2 := a2val;
              if rline[9] = '1' then btype := 'S';  // single
              if rline[9] = '2' then btype := 'D';  // double
              if rline[9] = '3' then btype := 'T';  // triple
              if rline[9] = '4' then btype := 'A';  // aromatic
              if rline[9] = '5' then btype := 'l';  // single or double
              if rline[9] = '6' then btype := 's';  // single or aromatic
              if rline[9] = '7' then btype := 'd';  // double or aromatic
              if rline[9] = '8' then btype := 'a';  // any
              if rline[9] = '9' then btype := 'a';  // any in JSME;  v0.5b
              if pos(btype,'lsda') > 0 then found_querymol := true;
              arom := false;
              q_arom := false;  // v0.3p
              // read aromaticity flag from CheckMol-tweaked MDL molfile
              if (btype = 'A') or (rline[8] = '0') then
                begin
                  arom := true;
                  if (rline[8] = '0') then found_arominfo := true;
                end;
              tmpstr := copy(rline,13,3);  // new in v0.3d: read ring_count from tweaked molfile
              val(tmpstr,rc,code);
              if ((code <> 0) or (rc < 0) or (progmode = pmCheckmol) or tmfmismatch) then 
                ring_count := 0 else ring_count := rc;  // v0.3n: added tmfmismatch check
              tmpstr := copy(rline,16,3);  // new in v0.3d: read bond topology;
              val(tmpstr,bt,code);         // extended features are encoded by leading zero
              if ((code <> 0) or (bt < 0) or (bt > 2)) then topo := btopo_any else   // v0.3n changed >5 into >2
                begin
                  if (tmpstr[2] = '0') then topo := bt + 3 else topo := bt;
                end;
              // new in v0.3d: add stereo property from MDL "stereo care" flag in atom block
              stereo := bstereo_any;
              if (btype ='D') then
                begin
                  if (atom^[a1].stereo_care and atom^[a2].stereo_care) then
                    begin                      // this is the MDL-conformant encoding,
                      stereo := bstereo_xyz;   // for an alternative see below
                      ez_flag := true;         // v0.3f
                    end else
                    begin
                      tmpstr := copy(rline,10,3);  // new in v0.3d: read bond stereo specification;
                      val(tmpstr,bs,code);         // this extended feature is encoded by a leading zero
                      mdl_stereo := bs;            // v0.3n
                      if ((code <> 0) or (bs <= 0) or (bs > 2)) then stereo := bstereo_any 
                        else stereo := bstereo_xyz;
                      if (tmpstr[2] = '0') then stereo := bstereo_xyz;
                    end;
                end;
              //if stereo <> bstereo_any then ez_search := true;
              if (stereo <> bstereo_any) then ez_flag := true;  // changed in v0.3f
              if (btype ='S') and (length(rline)>11) and (rline[12]='1') then stereo := bstereo_up;
              if (btype ='S') and (length(rline)>11) and (rline[12]='6') then stereo := bstereo_down;
              tmpstr := copy(rline,10,3);  // new in v0.3n: save original bond stereo specification;
              val(tmpstr,bs,code);         // v0.3n
              mdl_stereo := bs;            // v0.3n
            end;
          //if atom^[a1val].heavy and atom^[a2val].heavy then inc(n_heavybonds);  // moved down; v0.4b
        end;
    end;  // if (n_bonds > 0)...
  sepcount := 0;
  while (ri < molbufindex) and (sepcount < 1) do
    begin
      inc(ri);
      rline := molbuf^[ri];
      if (pos('M  CHG',rline) > 0) then 
        begin                   // new in v0.3j
          if clearcharges then  // "M  CHG" supersedes all "old-style" charge values
            begin
              for i := 1 to n_atoms do atom^[i].formal_charge := 0;
            end;
          read_charges(rline);
          clearcharges := false;  // subsequent "M  CHG" lines must not clear previous values
        end;
      if (pos('M  ISO',rline) > 0) then read_isotopes(rline);  // v0.3p
      if (pos('M  RAD',rline) > 0) then read_radicals(rline);  // v0.3p
      if (pos('$$$$',rline)>0) then
        begin
          inc(sepcount);
          sep_label := get_sep_label(rline);   // v0.4c
          if (molbufindex > (ri + 2)) then mol_in_queue := true;  // we assume this is an SDF file
        end;
    end;
  if (n_bonds > 0) then  // v0.4b  (must be done after "M  ISO" check)
    begin
      for n := 1 to n_bonds do
        begin
          a1val := bond^[n].a1;
          a2val := bond^[n].a2;
          if atom^[a1val].heavy and atom^[a2val].heavy then inc(n_heavybonds);
        end;
    end;    
  fillchar(ring^,sizeof(ringlist),0);
  for n := 1 to max_rings do  // new in v0.3
    begin
      ringprop^[n].size     := 0;
      ringprop^[n].arom     := false;
      ringprop^[n].envelope := false;
    end;
  li := ri + 1;
end;


procedure lblank(cols:integer; var nstr:string);
// inserts leading blanks up to a given number of total characters
begin
  if length(nstr) >= cols then exit;
  while length(nstr) < cols do nstr := ' '+nstr;
end;

function get_isotopemass(el:str2;nucnum:integer):double;  // v0.4d
var
  res : double;
begin
  res := nucnum;
  if (el = 'H ') and (nucnum = 1)  then res := 1.0078;
  if (el = 'H ') and (nucnum = 2)  then res := 2.0141;
  if (el = 'H ') and (nucnum = 3)  then res := 3.0160;
  if (el = 'C ') and (nucnum = 12) then res := 12.0000;
  if (el = 'C ') and (nucnum = 13) then res := 13.0034;
  if (el = 'C ') and (nucnum = 14) then res := 14.0032;
  if (el = 'F ') and (nucnum = 18) then res := 18.0009;
  if (el = 'F ') and (nucnum = 19) then res := 18.9984;
  if (el = 'Cl') and (nucnum = 35) then res := 34.9689;
  if (el = 'Cl') and (nucnum = 37) then res := 36.9659;
  if (el = 'Br') and (nucnum = 79) then res := 78.9183;
  if (el = 'Br') and (nucnum = 81) then res := 80.9163;
  if (el = 'I ') and (nucnum = 127) then res := 126.9045;
  if (el = 'I ') and (nucnum = 129) then res := 128.9050;
  if (el = 'N ') and (nucnum = 14) then res := 14.0031;
  if (el = 'N ') and (nucnum = 15) then res := 15.0001;
  if (el = 'O ') and (nucnum = 16) then res := 15.9949;
  if (el = 'O ') and (nucnum = 17) then res := 16.9991;
  if (el = 'S ') and (nucnum = 32) then res := 31.9721;
  if (el = 'S ') and (nucnum = 34) then res := 33.9679;
  if (el = 'P ') and (nucnum = 31) then res := 30.9738;
  // to be continued...
  get_isotopemass := res;
end;

procedure calc_mf_mw;   // v0.4d
const
  l1str : string = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
  l2str : string = ' abcdefghijklmnopqrstuvwxyz';
var
  i, j, k, l : integer;
  el, eltmp : str2;
  Htotal, Hexp, Himp, Himpsum, nucnum : integer;
  Hcount, Dcount, Tcount : integer;
  am, mwsum : double;
begin
  mol_formula := '';
  mol_weight  := 0;
  Hcount := 0;
  Dcount := 0;
  Tcount := 0;
  Himpsum := 0;
  mwsum := 0;
  for j := 1 to max_atomicnum do pt[j].count := 0;
  if (n_atoms > 0) then
    begin
      for i := 1 to n_atoms do
        begin
          el := atom^[i].element;
          el[2] := lowercase(el[2]);
          Himp := atom^[i].Htot - atom^[i].Hexp;
          if ((el = 'H ') or (el = 'D ') or (el = 'T ')) then
            begin
              if (el = 'D ') then inc(Dcount);
              if (el = 'T ') then inc(Tcount);
              if (el = 'H ') then
                begin
                  if (atom^[i].nucleon_number = 3) then inc(Tcount) else
                  if (atom^[i].nucleon_number = 2) then inc(Dcount) else inc(Hcount);
                end;            
            end else
            begin  // all other elements go here
              Himpsum := Himpsum + Himp;
              for j := 1 to max_atomicnum do
                begin
                  if (pt[j].el = el) then inc(pt[j].count);
                end;
            end;
        end;
      // now generate the molecular formula
      Hcount := Hcount + Himpsum;
      if (pt[6].count > 0) then mol_formula := 'C';
      if (pt[6].count > 1) then mol_formula := mol_formula + inttostr(pt[6].count);
      if (Hcount > 0) then mol_formula := mol_formula + 'H';
      if (Hcount > 1) then mol_formula := mol_formula + inttostr(Hcount);
      if (Dcount > 0) then mol_formula := mol_formula + 'D';
      if (Dcount > 1) then mol_formula := mol_formula + inttostr(Dcount);
      if (Tcount > 0) then mol_formula := mol_formula + 'T';
      if (Tcount > 1) then mol_formula := mol_formula + inttostr(Tcount);
      // and now the rest in alphabetical order
      eltmp := '  ';
      for k := 1 to length(l1str) do
        begin
          eltmp[1] := l1str[k];
          for l := 1 to length(l2str) do
            begin
              eltmp[2] := l2str[l];
              for j := 1 to max_atomicnum do
                begin
                  if ((eltmp <> 'C ') and (eltmp <> 'H ') and (eltmp <> 'T ')) then
                    begin
                      if ((eltmp = pt[j].el) and (pt[j].count > 0)) then
                        begin
                          if (eltmp[2] = ' ') then mol_formula := mol_formula + eltmp[1] else
                                                   mol_formula := mol_formula + eltmp;
                          if (pt[j].count > 1) then mol_formula := mol_formula + inttostr(pt[j].count);
                        end;
                    end;
                end;
            end;
        end;
      // finally, calculate the molecular weight
      for i := 1 to n_atoms do
        begin
          el := atom^[i].element;
          el[2] := lowercase(el[2]);
          if ((el <> 'H ') and (el <> 'D ') and (el <> 'T ')) then
            begin
              for j := 1 to max_atomicnum do
                begin
                  if (pt[j].el = el) then
                    begin
                      am := pt[j].am;
                      nucnum := atom^[i].nucleon_number;
                      if (nucnum > 0) then am := get_isotopemass(el,nucnum);
                      mwsum := mwsum + am;
                    end;
                end;
            end;
        end;
      mwsum := mwsum + Hcount*1.0079 + Dcount*get_isotopemass('H ',2) + Tcount*get_isotopemass('H ',3);
      mol_weight := mwsum;
    end;  // if n_atoms > 0 ...
end;

procedure write_MDLmolfile;
var
  i : integer;
  tmpstr : string;
  wline : string;
  a_chg, a_iso, a_rad : integer;
  tmflabel : string;  // v0.3m
  mfmwstr  : string;  // v0.4d
begin
  str(tweaklevel,tmflabel);                                   // v0.3m
  while (length(tmflabel) < 2) do tmflabel := '0' + tmflabel; // v0.3m
  tmflabel := 'TMF' + tmflabel;                               // v0.3m
  if length(molname)>80 then molname := copy(molname,1,80);
  writeln(molname);
  write('  CheckMol                        ',tmflabel);  // v0.3m
  if ringsearch_mode = rs_sar then write(':r0');         // v0.3m
  if ringsearch_mode = rs_ssr then write(':r1');         // v0.3m
  if opt_metalrings then write(':m1') else write(':m0'); // v0.3m
  writeln;
  str(mol_weight:1:2,mfmwstr);                            // v0.4d
  mfmwstr := 'MF=' + mol_formula + ' MW=' + mfmwstr ;  // v0.4d
  while (length(molcomment) > 0) and                   // v0.4d
        (molcomment[length(molcomment)] = ' ') do
    delete(molcomment,length(molcomment),1);
  write(molcomment);
  if (molcomment <> '') then write(' ');
  if (length(molcomment) + length(mfmwstr) < 80) and 
     (pos('MF=',molcomment) = 0) then write(mfmwstr);  // v0.4d
  writeln;   // v0.4d
  wline := ''; tmpstr := '';
  str(n_atoms:1,tmpstr); lblank(3,tmpstr);
  wline := wline+tmpstr; tmpstr := '';   // first 3 digits: number of atoms
  str(n_bonds:1,tmpstr); lblank(3,tmpstr);
  wline := wline+tmpstr; tmpstr := '';   // next 3 digits: number of bonds
  tmpstr := '  0';
  wline := wline + tmpstr; tmpstr := ''; // next 3 digits: number of atom lists (not used by us)
  {$IFDEF reduced_SAR}
  str(n_countablerings:1,tmpstr);  // v0.3n; changed n_rings into n_countablerings
  {$ELSE}
  str(n_rings:1,tmpstr);
  {$ENDIF}
  lblank(3,tmpstr);
  wline := wline+tmpstr; tmpstr := '';   // officially "obsolete", we use it for the number of rings
  wline := wline + '  ';                                            // v0.3n: obey chiral flag
  if chir_flag then wline := wline + '1' else wline := wline + '0'; // v0.3n
  wline := wline + '               999 V2000';                      // v0.3n (adjust string length)
  writeln(wline);
  for i := 1 to n_atoms do
    begin
      wline := '';
      str(atom^[i].X:1:4,tmpstr); lblank(10,tmpstr); wline := wline + tmpstr;
      str(atom^[i].Y:1:4,tmpstr); lblank(10,tmpstr); wline := wline + tmpstr;
      str(atom^[i].Z:1:4,tmpstr); lblank(10,tmpstr); wline := wline + tmpstr;
      tmpstr := atom^[i].element;
      tmpstr := lowercase(tmpstr);
      tmpstr[1] := upcase(tmpstr[1]);
      if (tmpstr = 'H ') and (atom^[i].nucleon_number = 2) and keep_DT then tmpstr := 'D '; // v0.3p
      if (tmpstr = 'H ') and (atom^[i].nucleon_number = 3) and keep_DT then tmpstr := 'T '; // v0.3p
      //wline := wline + ' '+atom^[i].element+' ';
      wline := wline + ' '+tmpstr+' ';
      wline := wline + ' 0';   // mass difference (isotopes)
      // now we code aromaticity into the old-style charge column (charges are now in the M  CHG line)
      if atom^[i].arom then tmpstr := ' 00' else tmpstr := '  0'; wline := wline + tmpstr;
      wline := wline + '  0  0  0  0  0  0  0  0  0  0';
      writeln(wline);
    end;
  for i := 1 to n_bonds do
    begin
      wline := '';
      str(bond^[i].a1:1,tmpstr); lblank(3,tmpstr); wline := wline + tmpstr;
      str(bond^[i].a2:1,tmpstr); lblank(3,tmpstr); wline := wline + tmpstr;
      if bond^[i].btype = 'S' then tmpstr := '  1';
      if bond^[i].btype = 'D' then tmpstr := '  2';
      if bond^[i].btype = 'T' then tmpstr := '  3';
      if bond^[i].btype = 'A' then tmpstr := '  4';
      if bond^[i].btype = 'l' then tmpstr := '  5';
      if bond^[i].btype = 's' then tmpstr := '  6';
      if bond^[i].btype = 'd' then tmpstr := '  7';
      if bond^[i].btype = 'a' then tmpstr := '  8';
      // now encode our own aromaticity information
      if bond^[i].arom then tmpstr[2] := '0';
      wline := wline + tmpstr;  // next, encode bond stereo property (v0.3f)
      //if (bond^[i].stereo = bstereo_up) then wline := wline + '  1' else
      //  if (bond^[i].stereo = bstereo_down) then wline := wline + '  6' else
      //    wline := wline + '  0';
      // restore original value from MDL molfile (v0.3n)
      wline := wline + '  ' + inttostr(bond^[i].mdl_stereo);    
      tmpstr := '';
      // now encode the ring_count of this bond (using a field which officially is "not used")
      tmpstr := inttostr(bond^[i].ring_count);
      while (length(tmpstr) < 3) do tmpstr := ' '+tmpstr;
      wline := wline + tmpstr+ '  0  0';
      writeln(wline);
    end;
  for i := 1 to n_atoms do
    begin
      a_chg := atom^[i].formal_charge;
      if a_chg <> 0 then
        begin
          wline := 'M  CHG  1 ';
          str(i:1,tmpstr); lblank(3,tmpstr); wline := wline + tmpstr +' ';
          str(a_chg:1,tmpstr); lblank(3,tmpstr); wline := wline + tmpstr;
          writeln(wline);
        end;
    end;
  for i := 1 to n_atoms do  // v0.3p
    begin
      a_iso := atom^[i].nucleon_number;
      if a_iso <> 0 then
        begin
          if (atom^[i].element <> 'H ') or (keep_DT = false) then
            begin
              wline := 'M  ISO  1 ';
              str(i:1,tmpstr); lblank(3,tmpstr); wline := wline + tmpstr +' ';
              str(a_iso:1,tmpstr); lblank(3,tmpstr); wline := wline + tmpstr;
              writeln(wline);
            end;
        end;
    end;
  for i := 1 to n_atoms do  // v0.3p
    begin
      a_rad := atom^[i].radical_type;
      if a_rad <> 0 then
        begin
          wline := 'M  RAD  1 ';
          str(i:1,tmpstr); lblank(3,tmpstr); wline := wline + tmpstr +' ';
          str(a_rad:1,tmpstr); lblank(3,tmpstr); wline := wline + tmpstr;
          writeln(wline);
        end;
    end;
  writeln('M  END');
end;


//============= chemical processing functions & procedures ============

procedure add2fgloc(row:integer;a:integer);  // v0.5
var
  i, n, p : integer;
  a_found : boolean;
begin
  if (row < 0) or (row > used_fg) then exit;
  n := fgloc^[row,0];
  a_found := false;
  if (n > 0) then
    begin
      for i := 1 to n do 
        begin
          if fgloc^[row,i] = a then a_found := true;
        end;
    end;
  if (a_found = false) and (n < max_fgpos) then
    begin
      inc(n);
      fgloc^[row,0] := n;
      fgloc^[row,n] := a;
    end;
end;


function is_electroneg(a_el:str2):boolean;
// new in v0.3j
var
  res : boolean;
begin
  res := false;;
  if a_el = 'N ' then res := true;
  if a_el = 'P ' then res := true;
  if a_el = 'O ' then res := true;
  if a_el = 'S ' then res := true;
  if a_el = 'SE' then res := true;
  if a_el = 'TE' then res := true;
  if a_el = 'F ' then res := true;
  if a_el = 'CL' then res := true;
  if a_el = 'BR' then res := true;
  if a_el = 'I ' then res := true;
  is_electroneg := res;
end;


procedure count_neighbors;
// counts heavy-atom neighbors and explicit hydrogens
var
  i : integer;
begin
  if (n_atoms < 1) or (n_bonds < 1) then exit;
  for i := 1 to n_bonds do
    begin
      if atom^[bond^[i].a1].heavy then inc(atom^[(bond^[i].a2)].neighbor_count);
      if atom^[bond^[i].a2].heavy then inc(atom^[(bond^[i].a1)].neighbor_count);
      if (atom^[(bond^[i].a1)].element = 'H ') then inc(atom^[(bond^[i].a2)].Hexp);
      if (atom^[(bond^[i].a2)].element = 'H ') then inc(atom^[(bond^[i].a1)].Hexp);
      // plausibility check (new in v02.i)
      if (atom^[(bond^[i].a1)].neighbor_count > max_neighbors) or 
         (atom^[(bond^[i].a2)].neighbor_count > max_neighbors) then
         begin
           mol_OK := false;
           //writeln('invalid molecule!');
         end;
    end;
end;


function get_neighbors(id:integer):neighbor_rec;
var
  i : integer;
  nb_tmp : neighbor_rec;
  nb_count : integer;
begin
  fillchar(nb_tmp,sizeof(neighbor_rec),0);
  nb_count := 0;
  for i := 1 to n_bonds do
    begin
      if ((bond^[i].a1 = id) and (nb_count < max_neighbors)) and (atom^[bond^[i].a2].heavy) then
        begin
          inc(nb_count);
          nb_tmp[nb_count] := bond^[i].a2;
        end;
      if ((bond^[i].a2 = id) and (nb_count < max_neighbors)) and (atom^[bond^[i].a1].heavy) then
        begin
          inc(nb_count);
          nb_tmp[nb_count] := bond^[i].a1;
        end;
    end;
  get_neighbors := nb_tmp;
end;


function get_ndl_neighbors(id:integer):neighbor_rec;
var
  i : integer;
  nb_tmp : neighbor_rec;
  nb_count : integer;  // v0.3i: use max_neighbors instead of a fixed value of 8
begin
  fillchar(nb_tmp,sizeof(neighbor_rec),0);
  nb_count := 0;
  for i := 1 to ndl_n_bonds do
    begin
      if ((ndl_bond^[i].a1 = id) and (nb_count < max_neighbors)) and (ndl_atom^[ndl_bond^[i].a2].heavy) then
        begin
          inc(nb_count);
          nb_tmp[nb_count] := ndl_bond^[i].a2;
        end;
      if ((ndl_bond^[i].a2 = id) and (nb_count < max_neighbors)) and (ndl_atom^[ndl_bond^[i].a1].heavy) then
        begin
          inc(nb_count);
          nb_tmp[nb_count] := ndl_bond^[i].a1;
        end;
    end;
  get_ndl_neighbors := nb_tmp;
end;


function get_nextneighbors(id:integer;prev_id:integer):neighbor_rec;
var
  i : integer;
  nb_tmp : neighbor_rec;
  nb_count : integer;
begin
  // gets all neighbors except prev_id (usually the atom where we came from
  fillchar(nb_tmp,sizeof(neighbor_rec),0);
  nb_count := 0;
  for i := 1 to n_bonds do
    begin
      if ((bond^[i].a1 = id) and (bond^[i].a2 <> prev_id) and (nb_count < max_neighbors)) 
        and (atom^[bond^[i].a2].heavy) then
        begin
          inc(nb_count);
          nb_tmp[nb_count] := bond^[i].a2;
        end;
      if ((bond^[i].a2 = id) and (bond^[i].a1 <> prev_id) and (nb_count < max_neighbors)) 
        and (atom^[bond^[i].a1].heavy) then
        begin
          inc(nb_count);
          nb_tmp[nb_count] := bond^[i].a1;
        end;
    end;
  get_nextneighbors := nb_tmp;
end;


function get_ndl_nextneighbors(id:integer;prev_id:integer):neighbor_rec;
var
  i : integer;
  nb_tmp : neighbor_rec;
  nb_count : integer;
begin
  // gets all neighbors except prev_id (usually the atom where we came from
  fillchar(nb_tmp,sizeof(neighbor_rec),0);
  nb_count := 0;
  for i := 1 to ndl_n_bonds do
    begin
      if ((ndl_bond^[i].a1 = id) and (ndl_bond^[i].a2 <> prev_id) and 
        (nb_count < max_neighbors)) and (ndl_atom^[ndl_bond^[i].a2].heavy) then
        begin
          inc(nb_count);
          nb_tmp[nb_count] := ndl_bond^[i].a2;
        end;
      if ((ndl_bond^[i].a2 = id) and (ndl_bond^[i].a1 <> prev_id) and 
        (nb_count < max_neighbors)) and (ndl_atom^[ndl_bond^[i].a1].heavy) then
        begin
          inc(nb_count);
          nb_tmp[nb_count] := ndl_bond^[i].a1;
        end;
    end;
  get_ndl_nextneighbors := nb_tmp;
end;


function path_pos(id:integer;a_path:ringpath_type):integer;  // new version in v0.3l
var
  i, pp : integer;
begin
  pp := 0;
  for i := 1 to max_ringsize do
    begin
      if (a_path[i] = id) then 
        begin
          pp := i;
          break;
        end;
    end;
  path_pos := pp;
end;


function path_length(a_path:ringpath_type):integer;
begin
  if (a_path[max_ringsize] <> 0) and (path_pos(0,a_path)=0) then path_length := max_ringsize else
    begin
      path_length := path_pos(0,a_path)-1;
    end;
end;


function matchpath_pos(id:integer;a_path:matchpath_type):integer;
var
  i, pp : integer;
begin
  pp := 0;
  for i := max_matchpath_length downto 1 do
    begin
      if (a_path[i] = id) then pp := i;
    end;
  matchpath_pos := pp;
end;


function matchpath_length(a_path:matchpath_type):integer;
begin
  if a_path[max_matchpath_length] <> 0 then matchpath_length := max_matchpath_length else
    begin
      matchpath_length := matchpath_pos(0,a_path)-1;
    end;
end;


function get_bond(ba1,ba2:integer):integer;
var
  i, b_id : integer;
begin
  b_id := 0;
  if n_bonds > 0 then begin
    for i := 1 to n_bonds do
      begin
        if ((bond^[i].a1 = ba1) and (bond^[i].a2 = ba2)) or
           ((bond^[i].a1 = ba2) and (bond^[i].a2 = ba1)) then
           b_id := i;
      end;
  end;
  get_bond := b_id;
end;


function get_ndl_bond(ba1,ba2:integer):integer;
var
  i, b_id : integer;
begin
  b_id := 0;
  if ndl_n_bonds > 0 then begin
    for i := 1 to ndl_n_bonds do
      begin
        if ((ndl_bond^[i].a1 = ba1) and (ndl_bond^[i].a2 = ba2)) or
           ((ndl_bond^[i].a1 = ba2) and (ndl_bond^[i].a2 = ba1)) then
           b_id := i;
      end;
  end;
  get_ndl_bond := b_id;
end;


procedure order_ringpath(var r_path:ringpath_type);
// order should be: array starts with atom of lowest number, followed by neighbor atom with lower number
var
  i, pl : integer;
  a_ref, a_left, a_right, a_tmp : integer;
begin
  pl := path_length(r_path);
  if (pl < 3) then exit;
  a_ref := n_atoms;  // start with highest possible value for an atom number
  for i := 1 to pl do
    begin
      if r_path[i] < a_ref then a_ref := r_path[i];  // find the minimum value ==> reference atom
    end;
  if a_ref < 1 then exit;  // just to be sure
  if path_pos(a_ref,r_path) < pl then a_right := r_path[(path_pos(a_ref,r_path)+1)] else a_right := r_path[1];
  if path_pos(a_ref,r_path) > 1 then a_left := r_path[(path_pos(a_ref,r_path)-1)] else a_left := r_path[pl];
  if a_right = a_left then exit;  // should never happen
  if a_right < a_left then
    begin  // correct ring numbering direction, only shift of the reference atom to the left end required
      while path_pos(a_ref,r_path) > 1 do
        begin
          a_tmp := r_path[1];
          for i := 1 to (pl - 1) do r_path[i] := r_path[(i+1)];
          r_path[pl] := a_tmp;
        end;
    end else
    begin  // wrong ring numbering direction, two steps required
      while path_pos(a_ref,r_path) < pl do
        begin  // step one: create "mirrored" ring path with reference atom at right end
          a_tmp := r_path[pl];
          for i := pl downto 2 do r_path[i] := r_path[(i-1)];
          r_path[1] := a_tmp;
        end;
      for i := 1 to (pl div 2) do
        begin  // one more mirroring
          a_tmp := r_path[i];
          r_path[i] := r_path[(pl+1)-i];
          r_path[(pl+1)-i] := a_tmp;
        end;
    end;
end;


function ringcompare(rp1,rp2:ringpath_type):integer;
var
  i, j, rc, rs1, rs2 : integer;
  n_common, max_cra : integer;
begin
  rc := 0;
  n_common := 0;
  rs1 := path_length(rp1);
  rs2 := path_length(rp2);
  if rs1 < rs2 then max_cra := rs1 else max_cra := rs2;
  for i := 1 to rs1 do
    for j := 1 to rs2 do
      if rp1[i] = rp2[j] then inc(n_common);
  if (rs1 = rs2) and (n_common = max_cra) then rc := 0 else
    begin
      if n_common = 0 then inc(rc,8);
      if n_common < max_cra then inc(rc,4) else
        begin
          if rs1 < rs2 then inc(rc,1) else inc(rc,2);
        end;
    end;
  ringcompare := rc;
end;


function rc_identical(rc_int:integer):boolean;
begin
  if rc_int = 0 then rc_identical := true else rc_identical := false;
end;


function rc_1in2(rc_int:integer):boolean;
begin
  if odd(rc_int) then rc_1in2 := true else rc_1in2 := false;
end;


function rc_2in1(rc_int:integer):boolean;
begin
  rc_int := rc_int div 2;
  if odd(rc_int) then rc_2in1 := true else rc_2in1 := false;
end;


function rc_different(rc_int:integer):boolean;
begin
  rc_int := rc_int div 4;
  if odd(rc_int) then rc_different := true else rc_different := false;
end;


function rc_independent(rc_int:integer):boolean;
begin
  rc_int := rc_int div 8;
  if odd(rc_int) then rc_independent := true else rc_independent := false;
end;


function is_newring(n_path:ringpath_type):boolean;
var
  i, j : integer;
  nr, same_ring : boolean;
  tmp_path : ringpath_type;
  rc_result : integer;
  found_ring : boolean;
  pl : integer;  // new in v0.3
begin
  nr := true;
  pl := path_length(n_path);  // new in v0.3
  if n_rings > 0 then
    begin
      case ringsearch_mode of
        rs_sar  : begin
                    found_ring := false;
                    i := 0;
                    while ((i < n_rings) and (not found_ring)) do
                      begin
                        inc(i);
                        if (pl = ringprop^[i].size) then  // compare only rings of same size
                          begin  
                            same_ring := true;   
                            for j := 1 to max_ringsize do
                              begin
                                if (ring^[i,j] <> n_path[j]) then same_ring := false;
                              end;
                            if same_ring then 
                              begin
                                nr := false;
                                found_ring := true;
                              end;
                          end;
                      end;  // while
                  end;
        rs_ssr  : begin
                    for i := 1 to n_rings do
                      begin
                        for j := 1 to max_ringsize do tmp_path[j] := ring^[i,j];
                        rc_result := ringcompare(n_path,tmp_path);
                        if rc_identical(rc_result) then nr := false;
                        if rc_1in2(rc_result) then
                          begin
                            // exchange existing ring by smaller one
                            for j := 1 to max_ringsize do ring^[i,j] := n_path[j];
                            // update ring property record (new in v0.3)
                            ringprop^[i].size := pl;
                            nr := false;
                            {$IFDEF debug}
                            debugoutput('replacing ring '+inttostr(i)+' by smaller one (ringsize: '+inttostr(path_length(n_path))+')');
                            {$ENDIF}
                          end;
                        if rc_2in1(rc_result) then
                          begin
                            // new ring contains existing one, but is larger ==> discard
                            nr := false;
                          end;
                      end;
                  end;
      end;  // case
    end;
  is_newring := nr;
end;


procedure add_ring(n_path:ringpath_type);
// new in v0.3: store rings in an ordered way (with ascending ring size)
var
  i, j, k, s, pl : integer;
  {$IFDEF debug}
  dstr : string;
  {$ENDIF}
begin
  pl := path_length(n_path);
  if pl < 1 then exit;
  if n_rings < max_rings then
    begin
      inc(n_rings);
      {$IFDEF debug}
      dstr := '';
      for j := 1 to pl do dstr := dstr + inttostr((n_path[j])) + ' ';
      debugoutput('adding ring '+inttostr(n_rings)+':  '+dstr);
      {$ENDIF}
      j := 0;
      if (n_rings > 1) then
        begin
          for i := 1 to (n_rings - 1) do
            begin
              s := ringprop^[i].size;
              if (pl >= s) then j := i;
            end;
        end;
      inc(j);  // the next position is ours
      if (j < n_rings) then
        begin  // push up the remaining rings by one position
          for k := n_rings downto (j+1) do
            begin
              ringprop^[k].size := ringprop^[(k-1)].size;
              //ringprop^[k].arom := ringprop^[(k-1)].arom;
              for i := 1 to max_ringsize do
                begin
                  ring^[k,i] := ring^[(k-1),i];
                end; 
            end;
        end;
      ringprop^[j].size := pl;  
      for i := 1 to max_ringsize do
        begin
          ring^[j,i] := n_path[i];
          //inc(atom^[(n_path[i])].ring_count);
        end; 
    end else 
    begin
      {$IFDEF debug}
      debugoutput('max_rings exceeded!');
      {$ENDIF}
    end;  
end;


function is_ringpath(s_path:ringpath_type):boolean;
var
  i, j : integer;
  nb : neighbor_rec;
  rp, new_atom : boolean;
  a_last, pl : integer;
  l_path : ringpath_type;
begin
  rp := false;
  new_atom := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  fillchar(l_path,sizeof(ringpath_type),0);
  pl := path_length(s_path);
  if pl < 1 then 
    begin 
      {$IFDEF debug}
      debugoutput('Oops! Got zero-length s_path!'); 
      {$ENDIF}
      exit; 
    end;
  for i := 1 to pl do
    begin
      l_path[i] := s_path[i];
    end;
  // check if the last atom is a metal and stop if opt_metalrings is not set (v0.3)
  if (opt_metalrings = false) then
    begin
      if atom^[l_path[pl]].metal then
        begin
          {$IFDEF debug}
          debugoutput('skipping metal in ring search'); 
          {$ENDIF}
          is_ringpath := false;
          exit;
        end;
    end;
  // check if ring is already closed
  if (pl > 2) and (l_path[pl] = l_path[1]) then
    begin
      l_path[pl] := 0;  // remove last entry (redundant!)
      order_ringpath(l_path);
      if is_newring(l_path) then
        begin
          if (n_rings < max_rings) then add_ring(l_path) else
            begin
              {$IFDEF debug}
              debugoutput('maximum number of rings exceeded!');
              {$ENDIF}
              is_ringpath := false;
              exit;
            end;
        end;
      rp := true;
      is_ringpath := true;
      exit;
    end;
  // any other case: ring is not (yet) closed
  a_last := l_path[pl];
  nb := get_neighbors(a_last);
  if atom^[a_last].neighbor_count > 1 then
    begin
      if ((rp = false) and (n_rings < max_rings)) then   // added in v0.2: check if max_rings is reached
        begin  // if ring is not closed, continue searching
          for i := 1 to atom^[a_last].neighbor_count do
            begin
              new_atom := true;
              for j := 2 to pl do if nb[i] = l_path[j] then 
                begin      // v0.3k
                  new_atom := false;
                  break;   // v0.3k
                end;
              // added in v0.1a: check if max_rings not yet reached
              // added in v0.2:  limit ring size to max_vringsize instead of max_ringsize
              if (new_atom) and (pl < max_vringsize) and (n_rings < max_rings) then
                begin
                  l_path[(pl+1)] := nb[i];
                  if (pl < max_ringsize-1) then l_path[pl+2] := 0;  // just to be sure
                  inc(recursion_level);                             // v0.3p (begin)
                  if (recursion_level > max_recursion_depth) then
                    begin
                      n_rings := max_rings;
                      is_ringpath := false;
                      exit;
                    end;                                            // v0.3p (end)
                  if is_ringpath(l_path) then rp := true;
                end;
            end;
        end;
    end;
  is_ringpath := rp;
end;


function is_ringbond(b_id:integer):boolean;
var
  i : integer;
  ra1, ra2 : integer;
  nb : neighbor_rec;
  search_path : ringpath_type;
  rb : boolean;
begin
  rb := false;
  recursion_level := 0;  // v0.3p
  ra1 := bond^[b_id].a1;
  ra2 := bond^[b_id].a2;
  fillchar(nb,sizeof(neighbor_rec),0);
  fillchar(search_path,sizeof(ringpath_type),0);
  nb := get_neighbors(ra2);
  if (atom^[ra2].neighbor_count > 1) and (atom^[ra1].neighbor_count > 1) then
    begin
      search_path[1] := ra1;
      search_path[2] := ra2;
      for i := 1 to atom^[ra2].neighbor_count do
        begin
          if (nb[i] <> ra1) and (atom^[nb[i]].heavy) then
            begin
              search_path[3] := nb[i];
              if is_ringpath(search_path) then rb := true;
            end;
        end;
    end;
  is_ringbond := rb;
end;


procedure chk_ringbonds;
var
  i : integer;
  a1rc, a2rc : integer;
begin
  if n_bonds < 1 then exit;
  for i := 1 to n_bonds do
    begin
      a1rc := atom^[(bond^[i].a1)].ring_count;
      a2rc := atom^[(bond^[i].a2)].ring_count;
      if ((n_rings = 0) or ((a1rc < n_rings) and (a2rc < n_rings) )) then
        begin
          if is_ringbond(i) then
            begin
              //inc(bond^[i].ring_count);
            end;
        end;
    end;
end;

(* v0.3d: moved procedure update_ringcount a bit down *)


function raw_hetbond_count(a:integer):integer;   // new in v0.2j, ignores bond order
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  hbc : integer;
begin
  hbc := 0;
  if (a > 0) and (a <= n_atoms) then
    begin
      fillchar(nb,sizeof(neighbor_rec),0);
      nb := get_neighbors(a);
      if atom^[a].neighbor_count > 0 then
        begin
          for i := 1 to atom^[a].neighbor_count do
            begin
              nb_el := atom^[(nb[i])].element;
              if (nb_el <> 'C ') and (nb_el <> 'A ') and (nb_el <> 'H ') and  // added 'D ' in v0.3n
                 (nb_el <> 'D ') and (nb_el <> 'LP') and (nb_el <> 'DU') then inc(hbc);
            end;
        end;
    end;
  raw_hetbond_count := hbc;
end;


function hetbond_count(a:integer):integer;
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  bt  : char;   // v0.4
  hbc : single;
  hbc_int : integer;  // v0.4
begin
  hbc := 0.0;
  if (a > 0) and (a <= n_atoms) then
    begin
      fillchar(nb,sizeof(neighbor_rec),0);
      nb := get_neighbors(a);
      if atom^[a].neighbor_count > 0 then
        begin
          for i := 1 to atom^[a].neighbor_count do
            begin
              nb_el := atom^[(nb[i])].element;
              if (nb_el <> 'C ') and (nb_el <> 'A ') and (nb_el <> 'H ') and
                 (nb_el <> 'D ') and (nb_el <> 'LP') and (nb_el <> 'DU') then  // added 'D ' in v0.3n
                begin
                  bt := bond^[(get_bond(a,nb[i]))].btype;  // v0.4
                  if (bt = 'S') then hbc := hbc + 1;
                  if (bt = 'A') then hbc := hbc + 1.5;
                  if (bt = 'D') then hbc := hbc + 2;
                  if (bt = 'T') then hbc := hbc + 3;
                end;
            end;
        end;
    end;
  hbc_int := round(hbc);  // v0.4
  hetbond_count := round(hbc_int);  // v0.4
end;


function hetatom_count(a:integer):integer;
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  hac : integer;
begin
  hac := 0;
  if (a > 0) and (a <= n_atoms) then
    begin
      fillchar(nb,sizeof(neighbor_rec),0);
      nb := get_neighbors(a);
      if atom^[a].neighbor_count > 0 then
        begin
          for i := 1 to atom^[a].neighbor_count do
            begin
              nb_el := atom^[(nb[i])].element;
              if (nb_el <> 'C ') and (nb_el <> 'H ') and (nb_el <> 'D ') and
                 (nb_el <> 'LP') and (nb_el <> 'DU') then inc(hac);  // added 'D ' in v0.3n
            end;
        end;
    end;
  hetatom_count := hac;
end;


function ndl_hetbond_count(a:integer):integer;
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  bt  : char;   // v0.4
  hbc : single;
  hbc_int : integer;  // v0.4
begin
  hbc := 0;
  if (a > 0) and (a <= n_atoms) then
    begin
      fillchar(nb,sizeof(neighbor_rec),0);
      nb := get_ndl_neighbors(a);
      if ndl_atom^[a].neighbor_count > 0 then
        begin
          for i := 1 to ndl_atom^[a].neighbor_count do
            begin
              nb_el := ndl_atom^[(nb[i])].element;
              if (nb_el <> 'C ') and (nb_el <> 'H ') and (nb_el <> 'D ') and
                 (nb_el <> 'LP') and (nb_el <> 'DU') then   // added 'D ' in v0.3n
                begin
                  bt := bond^[(get_ndl_bond(a,nb[i]))].btype;  // v0.4
                  if (bt = 'S') then hbc := hbc + 1;
                  if (bt = 'A') then hbc := hbc + 1.5;
                  if (bt = 'D') then hbc := hbc + 2;
                  if (bt = 'T') then hbc := hbc + 3;
                end;
            end;
        end;
    end;
  hbc_int := round(hbc);  // v0.4
  ndl_hetbond_count := round(hbc_int);  // v0.4
end;


function ndl_hetatom_count(a:integer):integer;
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  hac : integer;
begin
  hac := 0;
  if (a > 0) and (a <= ndl_n_atoms) then
    begin
      fillchar(nb,sizeof(neighbor_rec),0);
      nb := get_ndl_neighbors(a);
      if ndl_atom^[a].neighbor_count > 0 then
        begin
          for i := 1 to ndl_atom^[a].neighbor_count do
            begin  // note: query atoms like 'A' should be present only in the needle
              nb_el := ndl_atom^[(nb[i])].element;
              if (nb_el <> 'C ') and (nb_el <> 'A ') and (nb_el <> 'H ') and    // added 'D ' in v0.3n
                 (nb_el <> 'D ') and (nb_el <> 'LP') and (nb_el <> 'DU') then inc(hac);
            end;
        end;
    end;
  ndl_hetatom_count := hac;
end;


function is_oxo_C(id:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (id < 1) or (id > n_atoms) then exit;
  nb := get_neighbors(id);
  if (atom^[id].element = 'C ') and (atom^[id].neighbor_count > 0) then
    begin
      for i := 1 to atom^[id].neighbor_count do
        begin
          if (bond^[get_bond(id,nb[i])].btype = 'D') and
             ((atom^[(nb[i])].element = 'O ') { or
              (atom^[(nb[i])].element = 'S ')  or
              (atom^[(nb[i])].element = 'SE') } ) then     // no N, amidines are different...
             r := true;
        end;
    end;
  is_oxo_C := r;
end;


function is_thioxo_C(id:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (id < 1) or (id > n_atoms) then exit;
  nb := get_neighbors(id);
  if (atom^[id].element = 'C ') and (atom^[id].neighbor_count > 0) then
    begin
      for i := 1 to atom^[id].neighbor_count do
        begin
          if (bond^[get_bond(id,nb[i])].btype = 'D') and
             ((atom^[(nb[i])].element = 'S ')  or
              (atom^[(nb[i])].element = 'SE')) then     // no N, amidines are different...
             r := true;
        end;
    end;
  is_thioxo_C := r;
end;


function is_imino_C(id:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (id < 1) or (id > n_atoms) then exit;
  nb := get_neighbors(id);
  if (atom^[id].element = 'C ') and (atom^[id].neighbor_count > 0) then
    begin
      for i := 1 to atom^[id].neighbor_count do
        begin
          if (bond^[get_bond(id,nb[i])].btype = 'D') and
             (atom^[(nb[i])].element = 'N ') then
             r := true;
        end;
    end;
  is_imino_C := r;
end;


function is_true_imino_C(id:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  nb_el : str2;
  a_n, b : integer;  // v0.3j
begin
  r := false;
  a_n := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (id < 1) or (id > n_atoms) then exit;
  nb := get_neighbors(id);
  if (atom^[id].element = 'C ') and (atom^[id].neighbor_count > 0) then
    begin
      for i := 1 to atom^[id].neighbor_count do
        begin
          b := get_bond(id,nb[i]);  // v0.3j
          if (bond^[b].btype = 'D') and (bond^[b].arom = false) and  // v0.3j
             (atom^[(nb[i])].element = 'N ')  then a_n := nb[i];
        end;
      if (a_n > 0) then
        begin
          r := true;
          fillchar(nb,sizeof(neighbor_rec),0);
          nb := get_neighbors(a_n);
          for i := 1 to atom^[a_n].neighbor_count do
            begin
              nb_el := atom^[(nb[i])].element;
              if (nb_el <> 'C ') and (nb_el <> 'H ') and (nb_el <> 'D ') then r := false; // v0.3n: D
            end;  
        end;
    end;
  is_true_imino_C := r;
end;

function is_true_exocyclic_imino_C(id,r_id:integer):boolean;  // v0.3j
var
  i,j  : integer;
  r    : boolean;
  nb   : neighbor_rec;
  testring : ringpath_type;
  ring_size : integer;
  b : integer;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (id < 1) or (id > n_atoms) then exit;
  nb := get_neighbors(id);
  fillchar(testring,sizeof(ringpath_type),0);
  ring_size := ringprop^[r_id].size;  // v0.3j
  for j := 1 to ring_size do testring[j] := ring^[r_id,j];  // v0.3j
  if (atom^[id].element = 'C ') and (atom^[id].neighbor_count > 0) then
    begin
      for i := 1 to atom^[id].neighbor_count do
        begin
          b := get_bond(id,nb[i]);
          if (bond^[b].btype = 'D') and (bond^[b].arom = false) and
             (atom^[(nb[i])].element = 'N ') then
               begin
                 r := true;
                 for j := 1 to ring_size do
                   if nb[i] = ring^[r_id,j] then r := false;
               end;
        end;
    end;
  is_true_exocyclic_imino_C := r;
end;


function is_exocyclic_imino_C(id,r_id:integer):boolean;
var
  i,j  : integer;
  r    : boolean;
  nb   : neighbor_rec;
  testring : ringpath_type;
  ring_size : integer;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (id < 1) or (id > n_atoms) then exit;
  nb := get_neighbors(id);
  fillchar(testring,sizeof(ringpath_type),0);
  ring_size := ringprop^[r_id].size;  // v0.3j
  for j := 1 to ring_size do testring[j] := ring^[r_id,j];  // v0.3j
  if (atom^[id].element = 'C ') and (atom^[id].neighbor_count > 0) then
    begin
      for i := 1 to atom^[id].neighbor_count do
        begin
          if (bond^[get_bond(id,nb[i])].btype = 'D') and
             (atom^[(nb[i])].element = 'N ') then
               begin
                 r := true;
                 for j := 1 to ring_size do
                   if nb[i] = ring^[r_id,j] then r := false;
               end;
        end;
    end;
  is_exocyclic_imino_C := r;
end;


function find_exocyclic_methylene_C(id,r_id:integer):integer; 
var                    // renamed and rewritten in v0.3j
  i,j  : integer;
  r    : integer;
  nb   : neighbor_rec;
  testring : ringpath_type;
  ring_size : integer;
begin
  r := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (id < 1) or (id > n_atoms) then 
    begin
      find_exocyclic_methylene_C := 0;
      exit;
    end;
  nb := get_neighbors(id);
  fillchar(testring,sizeof(ringpath_type),0);
  ring_size := ringprop^[r_id].size;  // v0.3j
  for j := 1 to ring_size do testring[j] := ring^[r_id,j];  // v0.3j
  if (atom^[id].element = 'C ') and (atom^[id].neighbor_count > 0) then
    begin
      for i := 1 to atom^[id].neighbor_count do
        begin
          if (bond^[get_bond(id,nb[i])].btype = 'D') and
             (atom^[(nb[i])].element = 'C ') then
               begin
                 r := nb[i];
                 for j := 1 to ring_size do
                   if nb[i] = ring^[r_id,j] then r := 0;
               end;
        end;
    end;
  find_exocyclic_methylene_C := r;
end;


function is_hydroxy(a_view,a_ref:integer):boolean;
var
  r  : boolean;
begin
  r := false;
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'O3 ') and (atom^[a_ref].neighbor_count = 1) then r := true;
    end;
  is_hydroxy := r;
end;


function is_sulfanyl(a_view,a_ref:integer):boolean;
var
  r  : boolean;
begin
  r := false;
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'S3 ') and (atom^[a_ref].neighbor_count = 1) then r := true;
    end;
  is_sulfanyl := r;
end;


function is_amino(a_view,a_ref:integer):boolean;
var
  r  : boolean;
begin
  r := false;
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if ((atom^[a_ref].atype = 'N3 ') or (atom^[a_ref].atype = 'N3+'))
          and (atom^[a_ref].neighbor_count = 1) then r := true;
    end;
  is_amino := r;
end;


function is_alkyl(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  nb_el : str2;
  het_count : integer;
begin
  r := false;
  het_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') and
     (atom^[a_ref].atype = 'C3 ') and (atom^[a_ref].arom = false) then
    begin
      nb := get_nextneighbors(a_ref,a_view);
      for i := 1 to (atom^[a_ref].neighbor_count - 1) do
        begin
          nb_el := atom^[(nb[i])].element;
          if (nb_el <> 'C ') and (nb_el <> 'H ') and (nb_el <> 'D ') and
             (nb_el <> 'DU') and (nb_el <> 'LP') then inc(het_count);  // added 'D ' in v0.3n
        end;
      if het_count <= 1 then r := true;  // we consider (e.g.) alkoxyalkyl groups as alkyl
    end;
  is_alkyl := r;
end;


function is_true_alkyl(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  nb_el : str2;
  het_count : integer;
begin
  r := false;
  het_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') and
     (atom^[a_ref].atype = 'C3 ') and (atom^[a_ref].arom = false) then
    begin
      nb := get_nextneighbors(a_ref,a_view);
      for i := 1 to (atom^[a_ref].neighbor_count - 1) do
        begin
          nb_el := atom^[(nb[i])].element;
          if (nb_el <> 'C ') and (nb_el <> 'H ') and (nb_el <> 'D ') and
             (nb_el <> 'DU') then inc(het_count);  // added 'D ' in v0.3n
        end;
      if het_count = 0 then r := true;  //
    end;
  is_true_alkyl := r;
end;

function is_alkenyl(a_view,a_ref:integer):boolean;  // new in v0.3j; rewritten in v0.4b
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  nb_el : str2;
  nb_bond : integer;       // v0.4b
  nb_bt : char;            // v0.4b
  cc_dbl_count : integer;  // v0.4b
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') and
     (atom^[a_ref].atype = 'C2 ') and (atom^[a_ref].arom = false) then
    begin
      nb := get_nextneighbors(a_ref,a_view);
      cc_dbl_count := 0;
      for i := 1 to (atom^[a_ref].neighbor_count - 1) do
        begin
          nb_el := atom^[(nb[i])].element;
          nb_bond := get_bond(a_ref,nb[i]);
          if (bond^[nb_bond].btype = 'D') and (nb_el = 'C ') then inc(cc_dbl_count);
        end;
      if cc_dbl_count > 0 then r := true;       // we consider (e.g.) alkoxyalkenyl groups as alkenyl
    end;                                        // v0.3k: changed c2_count = 1 into c2_count >= 1
  is_alkenyl := r;
end;

function is_alkynyl(a_view,a_ref:integer):boolean;  // new in v0.3j
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  nb_at : str3;
  c1_count  : integer;
begin
  r := false;
  c1_count  := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') and
     (atom^[a_ref].atype = 'C1 ') and (atom^[a_ref].arom = false) then
    begin
      nb := get_nextneighbors(a_ref,a_view);
      for i := 1 to (atom^[a_ref].neighbor_count - 1) do
        begin
          nb_at := atom^[(nb[i])].atype;
          if (nb_at = 'C1 ') then inc(c1_count);  
        end;
      if (c1_count = 1) then r := true;
    end;
  is_alkynyl := r;
end;


function is_aryl(a_view,a_ref:integer):boolean;
var
  r  : boolean;
begin
  r := false;
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') and
     (atom^[a_ref].element = 'C ') and (atom^[a_ref].arom = true) then r := true;
  is_aryl := r;
end;


function is_alkoxy(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'O3 ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_alkyl(a_ref,nb[1]) then r := true;
        end;
    end;
  is_alkoxy := r;
end;


function is_siloxy(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'O3 ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if (atom^[(nb[1])].element = 'SI') then r := true;
        end;
    end;
  is_siloxy := r;
end;


function is_true_alkoxy(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'O3 ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_true_alkyl(a_ref,nb[1]) then r := true;
        end;
    end;
  is_true_alkoxy := r;
end;


function is_aryloxy(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'O3 ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_aryl(a_ref,nb[1]) then r := true;
        end;
    end;
  is_aryloxy := r;
end;

function is_alkenyloxy(a_view,a_ref:integer):boolean;  // v0.3j
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'O3 ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_alkenyl(a_ref,nb[1]) then r := true;
        end;
    end;
  is_alkenyloxy := r;
end;

function is_alkynyloxy(a_view,a_ref:integer):boolean;  // v0.3j
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'O3 ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_alkynyl(a_ref,nb[1]) then r := true;
        end;
    end;
  is_alkynyloxy := r;
end;


function is_alkylsulfanyl(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'S3 ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_alkyl(a_ref,nb[1]) then r := true;
        end;
    end;
  is_alkylsulfanyl := r;
end;


function is_true_alkylsulfanyl(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'S3 ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_true_alkyl(a_ref,nb[1]) then r := true;
        end;
    end;
  is_true_alkylsulfanyl := r;
end;


function is_arylsulfanyl(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'S3 ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_aryl(a_ref,nb[1]) then r := true;
        end;
    end;
  is_arylsulfanyl := r;
end;

function is_alkenylsulfanyl(a_view,a_ref:integer):boolean;  // v0.3j
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'S3 ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_alkenyl(a_ref,nb[1]) then r := true;
        end;
    end;
  is_alkenylsulfanyl := r;
end;

function is_alkynylsulfanyl(a_view,a_ref:integer):boolean;  // v0.3j
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].atype = 'S3 ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_alkynyl(a_ref,nb[1]) then r := true;
        end;
    end;
  is_alkynylsulfanyl := r;
end;

function is_alkylamino(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
  alkyl_count : integer;
begin
  r := false;
  alkyl_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_alkyl(a_ref,nb[1]) then inc(alkyl_count);
          if alkyl_count = 1 then  r := true  
        end;
    end;
  is_alkylamino := r;
end;


function is_dialkylamino(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  alkyl_count : integer;
begin
  r := false;
  alkyl_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if is_alkyl(a_ref,nb[i]) then inc(alkyl_count);
            end;
          if alkyl_count = 2 then  r := true  
        end;
    end;
  is_dialkylamino := r;
end;

function is_arylamino(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
  aryl_count : integer;
begin
  r := false;
  aryl_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_aryl(a_ref,nb[1]) then inc(aryl_count);
          if aryl_count = 1 then  r := true  
        end;
    end;
  is_arylamino := r;
end;


function is_diarylamino(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  aryl_count : integer;
begin
  r := false;
  aryl_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if is_aryl(a_ref,nb[i]) then inc(aryl_count);
            end;
          if aryl_count = 2 then  r := true
        end;
    end;
  is_diarylamino := r;
end;


function is_alkylarylamino(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  alkyl_count : integer;
  aryl_count  : integer;
begin
  r := false;
  alkyl_count := 0;
  aryl_count  := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if is_alkyl(a_ref,nb[i]) then inc(alkyl_count);
              if is_aryl(a_ref,nb[i]) then inc(aryl_count);
            end;
          if (alkyl_count = 1) and (aryl_count = 1) then  r := true  
        end;
    end;
  is_alkylarylamino := r;
end;

function is_C_monosubst_amino(a_view,a_ref:integer):boolean;  // new in v0.3j
// a_ref = N
var
  r  : boolean;
  nb : neighbor_rec;
  c_count : integer;
begin
  r := false;
  c_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if ((atom^[a_ref].atype = 'N3 ') or (atom^[a_ref].atype = 'NAM')) and 
        (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if (atom^[(nb[1])].element = 'C ') then inc(c_count);
          if (c_count = 1) then  r := true  
        end;
    end;
  is_C_monosubst_amino := r;
end;

function is_C_disubst_amino(a_view,a_ref:integer):boolean;  // new in v0.3j
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  b  : integer;
  c_count : integer;
begin
  r := false;
  c_count := 0;
  b := get_bond(a_view,a_ref);
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[b].btype = 'S') and (bond^[b].arom = false) then
    begin
      if ((atom^[a_ref].atype = 'N3 ') or (atom^[a_ref].atype = 'NAM')) and 
        (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if (atom^[(nb[i])].element = 'C ') then inc(c_count);
            end;   
          if (c_count = 2) then  r := true  
        end;
    end;
  is_C_disubst_amino := r;
end;


function is_subst_amino(a_view,a_ref:integer):boolean;
var
  r : boolean;
begin
  r := false;
  if (is_amino(a_view,a_ref)) or (is_alkylamino(a_view,a_ref)) or
     (is_arylamino(a_view,a_ref)) or (is_dialkylamino(a_view,a_ref)) or
     (is_alkylarylamino(a_view,a_ref)) or (is_diarylamino(a_view,a_ref)) then r := true;
  is_subst_amino := r;
end;


function is_true_alkylamino(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
  alkyl_count : integer;
begin
  r := false;
  alkyl_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if ((atom^[a_ref].atype = 'N3 ') or (atom^[a_ref].atype = 'N3+'))
         and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_true_alkyl(a_ref,nb[1]) then inc(alkyl_count);
          if alkyl_count = 1 then  r := true
        end;
    end;
  is_true_alkylamino := r;
end;


function is_true_dialkylamino(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  alkyl_count : integer;
begin
  r := false;
  alkyl_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if ((atom^[a_ref].atype = 'N3 ') or (atom^[a_ref].atype = 'N3+'))
         and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if is_true_alkyl(a_ref,nb[i]) then inc(alkyl_count);
            end;
          if alkyl_count = 2 then  r := true  
        end;
    end;
  is_true_dialkylamino := r;
end;


function is_true_alkylarylamino(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  alkyl_count : integer;
  aryl_count  : integer;
begin
  r := false;
  alkyl_count := 0;
  aryl_count  := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if ((atom^[a_ref].atype = 'N3 ') or (atom^[a_ref].atype = 'N3+'))
         and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if is_true_alkyl(a_ref,nb[i]) then inc(alkyl_count);
              if is_aryl(a_ref,nb[i]) then inc(aryl_count);
            end;
          if (alkyl_count = 1) and (aryl_count = 1) then  r := true  
        end;
    end;
  is_true_alkylarylamino := r;
end;


function is_hydroxylamino(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  oh_count : integer;
  het_count : integer;  // v0.3k
  nb_el : str2;         // v0.3k
begin
  r := false;
  oh_count  := 0;
  het_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count >= 2) then  // v0.3c
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to (atom^[a_ref].neighbor_count - 1) do  // v0.3c
            begin
              if is_hydroxy(a_ref,nb[i]) then inc(oh_count);
              nb_el := atom^[(nb[i])].element;                // v0.3k
              if (nb_el <> 'C ') and (nb_el <> 'H ') and      // v0.3k
                 (nb_el <> 'D ') and (nb_el <> 'DU') and (nb_el <> 'LP') then inc(het_count); // v0.3n: D
            end;
          if (oh_count = 1) and (het_count = 1) then  r := true  
        end;
    end;
  is_hydroxylamino := r;
end;


function is_nitro(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  o_count : integer;
  bond_count : integer;
begin
  r := false;
  o_count := 0;
  bond_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if (atom^[(nb[i])].element = 'O ') then inc(o_count);
              if (bond^[get_bond(a_ref,nb[i])].btype = 'S') then inc(bond_count);
              if (bond^[get_bond(a_ref,nb[i])].btype = 'D') then inc(bond_count,2);
            end;
          if (o_count = 2) and (bond_count >= 3) then  r := true  
        end;
    end;
  is_nitro := r;
end;


function is_azido(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
  bond_count : integer;
  n1, n2, n3 : integer;
begin
  r := false;
  bond_count := 0;
  n1 := 0; n2 := 0; n3 := 0;
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          n1 := a_ref;
          fillchar(nb,sizeof(neighbor_rec),0);
          nb := get_nextneighbors(n1,a_view);
          if (atom^[(nb[1])].element = 'N ') then
            begin
              n2 := nb[1];
              if (bond^[get_bond(n1,n2)].btype = 'S') then inc(bond_count);
              if (bond^[get_bond(n1,n2)].btype = 'D') then inc(bond_count,2);
              if (bond^[get_bond(n1,n2)].btype = 'T') then inc(bond_count,3);
            end;
          if (n2 > 0) and (atom^[n2].neighbor_count = 2) then
            begin
              fillchar(nb,sizeof(neighbor_rec),0);
              nb := get_nextneighbors(n2,n1);
              if (atom^[(nb[1])].element = 'N ') then
                begin
                  n3 := nb[1];
                  if (bond^[get_bond(n2,n3)].btype = 'S') then inc(bond_count);
                  if (bond^[get_bond(n2,n3)].btype = 'D') then inc(bond_count,2);
                  if (bond^[get_bond(n2,n3)].btype = 'T') then inc(bond_count,3);
                end;
            end;  
          if (n1 > 0) and (n2 > 0) and (n3 > 0) and (atom^[n3].neighbor_count = 1) and
             (bond_count > 3) then r := true  
        end;
    end;
  is_azido := r;
end;


function is_diazonium(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
  bond_count : integer;
  chg_count : integer;
  n1, n2 : integer;
begin
  r := false;
  bond_count := 0;
  chg_count := 0;
  n1 := 0; n2 := 0;
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          n1 := a_ref;
          chg_count := atom^[n1].formal_charge;
          fillchar(nb,sizeof(neighbor_rec),0);
          nb := get_nextneighbors(n1,a_view);
          if (atom^[(nb[1])].element = 'N ') then
            begin
              n2 := nb[1];
              chg_count := chg_count + atom^[n2].formal_charge;                      
              if (bond^[get_bond(n1,n2)].btype = 'S') then inc(bond_count);
              if (bond^[get_bond(n1,n2)].btype = 'D') then inc(bond_count,2);
              if (bond^[get_bond(n1,n2)].btype = 'T') then inc(bond_count,3);
            end;
          if (n1 > 0) and (n2 > 0) and (atom^[n2].neighbor_count = 1) and
             (bond_count >= 2) and (chg_count > 0) then r := true
        end;
    end;
  is_diazonium := r;
end;


function is_hydroximino_C(id:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  a_het : integer;
begin
  r := false;
  a_het := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (id < 1) or (id > n_atoms) then exit;
  nb := get_neighbors(id);
  if (atom^[id].element = 'C ') and (atom^[id].neighbor_count > 0) then
    begin
      for i := 1 to atom^[id].neighbor_count do
        begin
          if (bond^[get_bond(id,nb[i])].btype = 'D') and
             (atom^[(nb[i])].element = 'N ') and
             (hetbond_count(nb[i]) = 3) then
             a_het := nb[i];;
        end;
      if (a_het > 0) then
        begin
          fillchar(nb,sizeof(neighbor_rec),0);
          nb := get_neighbors(a_het);
          if (atom^[a_het].element = 'N ') and (atom^[a_het].neighbor_count > 0) then
            begin
              for i := 1 to atom^[a_het].neighbor_count do
                begin
                  if is_hydroxy(a_het,nb[i]) then r := true;
                end;
            end;
        end;
    end;
  is_hydroximino_C := r;
end;


function is_hydrazono_C(id:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  a_het : integer;
begin
  r := false;
  a_het := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (id < 1) or (id > n_atoms) then exit;
  nb := get_neighbors(id);
  if (atom^[id].element = 'C ') and (atom^[id].neighbor_count > 0) then
    begin
      for i := 1 to atom^[id].neighbor_count do
        begin
          if (bond^[get_bond(id,nb[i])].btype = 'D') and
             (atom^[(nb[i])].element = 'N ') { and
             (hetbond_count(nb[i]) = 3)  } then
             a_het := nb[i];;
        end;
      if (a_het > 0) then
        begin
          fillchar(nb,sizeof(neighbor_rec),0);
          nb := get_neighbors(a_het);
          if (atom^[a_het].element = 'N ') and (atom^[a_het].neighbor_count > 0) then
            begin
              for i := 1 to atom^[a_het].neighbor_count do
                begin
                  if (is_amino(a_het,nb[i])) or
                     (is_alkylamino(a_het,nb[i])) or
                     (is_alkylarylamino(a_het,nb[i])) or
                     (is_arylamino(a_het,nb[i])) or
                     (is_dialkylamino(a_het,nb[i])) or
                     (is_diarylamino(a_het,nb[i]))
                   then r := true;
                end;
            end;
        end;
    end;
  is_hydrazono_C := r;
end;


function is_alkoxycarbonyl(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (is_oxo_c(a_ref)) and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if is_alkoxy(a_ref,nb[i]) then r := true;
            end;
        end;    
    end;
  is_alkoxycarbonyl := r;  
end;


function is_aryloxycarbonyl(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (is_oxo_c(a_ref)) and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if is_aryloxy(a_ref,nb[i]) then r := true;
            end;
        end;
    end;
  is_aryloxycarbonyl := r;
end;


function is_carbamoyl(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (is_oxo_c(a_ref)) and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if (atom^[(nb[i])].atype = 'N3 ') or 
                 (atom^[(nb[i])].atype = 'NAM') then r := true;
            end;
        end;    
    end;
  is_carbamoyl := r;  
end;


function is_alkoxythiocarbonyl(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (is_thioxo_c(a_ref)) and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if is_alkoxy(a_ref,nb[i]) then r := true;
            end;
        end;    
    end;
  is_alkoxythiocarbonyl := r;  
end;


function is_aryloxythiocarbonyl(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (is_thioxo_c(a_ref)) and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if is_aryloxy(a_ref,nb[i]) then r := true;
            end;
        end;    
    end;
  is_aryloxythiocarbonyl := r;  
end;


function is_thiocarbamoyl(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (is_thioxo_c(a_ref)) and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if (atom^[(nb[i])].atype = 'N3 ') or
                 (atom^[(nb[i])].atype = 'NAM') then r := true;
            end;
        end;
    end;
  is_thiocarbamoyl := r;
end;


function is_alkanoyl(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (is_oxo_c(a_ref)) and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if (is_alkyl(a_ref,nb[i])) then r := true;
            end;
        end;    
    end;
  is_alkanoyl := r;  
end;


function is_aroyl(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (is_oxo_c(a_ref)) and (atom^[a_ref].neighbor_count = 3) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to 2 do
            begin
              if (is_aryl(a_ref,nb[i])) then r := true;
            end;
        end;    
    end;
  is_aroyl := r;  
end;


function is_acyl(a_view,a_ref:integer):boolean;
var
  r  : boolean;
begin
  r := false;
  if (is_alkanoyl(a_view,a_ref)) or (is_aroyl(a_view,a_ref)) then r := true;
  is_acyl := r;
end;

function is_acyl_gen(a_view,a_ref:integer):boolean;  // new in v0.3j
var
  r  : boolean;
begin
  r := false;
  if (is_oxo_C(a_ref)) then r := true;
  is_acyl_gen := r;
end;


function is_acylamino(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
  acyl_count : integer;
begin
  r := false;
  acyl_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          if is_acyl(a_ref,nb[1]) then inc(acyl_count);
          if acyl_count = 1 then  r := true  
        end;
    end;
  is_acylamino := r;
end;

function is_subst_acylamino(a_view,a_ref:integer):boolean;
var  // may be substituted _or_ unsubstituted acylamino group!
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  acyl_count : integer;
begin
  r := false;
  acyl_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count >= 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to (atom^[a_ref].neighbor_count - 1) do
            begin
              if is_acyl_gen(a_ref,nb[i]) then inc(acyl_count);  // v0.3j
            end;
          if (acyl_count) > 0 then r := true  
        end;
    end;
  is_subst_acylamino := r;
end;


function is_hydrazino(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  nr_count : integer;
begin
  r := false;
  nr_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count >= 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to (atom^[a_ref].neighbor_count - 1) do  // fixed in v0.3c
            begin
              if (is_amino(a_ref,nb[i])) or
                 (is_subst_amino(a_ref,nb[i])) then inc(nr_count);
            end;
          if (nr_count = 1) then r := true  
        end;
    end;
  is_hydrazino := r;
end;

procedure fgloc_set_hydrazino(a_ref:integer);   // v0.5
var
  i  : integer;
  nb : neighbor_rec;
  a_n : integer;
  b_id : integer;
begin
  if not opt_pos then exit;
  a_n := 0;  // v0.5
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_ref].heavy) then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count >= 2) then
        begin
          nb := get_neighbors(a_ref);
          for i := 1 to (atom^[a_ref].neighbor_count) do
            begin
              if (is_amino(a_ref,nb[i])) or
                 (is_subst_amino(a_ref,nb[i])) then a_n := nb[i];
            end;
          if a_n > 0 then
            begin
              b_id := get_bond(a_ref,a_n);
              add2fgloc(fg_hydrazine,b_id);
            end;
        end;
    end;
end;


function is_nitroso(a_view,a_ref:integer):boolean;  // new in v0.3j
var
  r  : boolean;
  nb : neighbor_rec;
  o_count : integer;
  a2 : integer;
begin
  r := false;
  o_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count = 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          a2 := nb[1];
          if (atom^[a2].element = 'O ') and (bond^[get_bond(a_ref,a2)].btype = 'D') then inc(o_count);
          if (o_count = 1) then r := true  
        end;
    end;
  is_nitroso := r;
end;


function is_subst_hydrazino(a_view,a_ref:integer):boolean;
var
  i  : integer;
  r  : boolean;
  nb : neighbor_rec;
  nr_count : integer;
  a2 : integer;
begin
  r := false;
  nr_count := 0;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_view].heavy) and (bond^[get_bond(a_view,a_ref)].btype = 'S') then
    begin
      if (atom^[a_ref].element = 'N ') and (atom^[a_ref].neighbor_count >= 2) then
        begin
          nb := get_nextneighbors(a_ref,a_view);
          for i := 1 to (atom^[a_ref].neighbor_count - 1) do
            begin
              a2 := nb[i];
              if (atom^[a2].element = 'N ') and
                 (not is_nitroso(a_ref,a2)) then inc(nr_count);  // v0.3j
            end;
          if (nr_count = 1) then r := true  
        end;
    end;
  is_subst_hydrazino := r;
end;


function is_cyano(a_view,a_ref:integer):boolean;
// a_view = C, a_ref = N
var
  r  : boolean;
begin
  r := false;
  if (atom^[a_view].atype = 'C1 ') and (bond^[get_bond(a_view,a_ref)].btype = 'T') and
     (atom^[a_ref].atype = 'N1 ') and (atom^[a_ref].neighbor_count = 1) then r := true;
  is_cyano := r;
end;


function is_cyano_c(a_ref:integer):boolean;
var
  i : integer;
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_ref].atype = 'C1 ') and (atom^[a_ref].neighbor_count > 0) then
    begin
      nb := get_neighbors(a_ref);
      for i := 1 to atom^[a_ref].neighbor_count do
        begin
          if (is_cyano(a_ref,nb[i])) then r := true;
        end;
    end;
  is_cyano_c := r;
end;


function is_nitrile(a_view,a_ref:integer):boolean;
// a_view = C, a_ref = N
var
  r  : boolean;
  nb : neighbor_rec;
  nb_el : str2;
begin
  r := false;
  if is_cyano(a_view,a_ref) then
  begin
    if (atom^[a_view].neighbor_count = 1) and
       (atom^[a_view].formal_charge = 0) then r := true else   // HCN is also a nitrile!
      begin
        nb := get_nextneighbors(a_view,a_ref);
        nb_el := atom^[(nb[1])].element;
        if (nb_el = 'C ') or (nb_el = 'H ') or (nb_el = 'D ') then r := true;  // v0.3n: D
      end;
    end;
  is_nitrile := r;
end;


function is_isonitrile(a_view,a_ref:integer):boolean;   // only recognized with CN triple bond!
// a_view = C, a_ref = N
var
  r  : boolean;
begin
  r := false;
  if (atom^[a_view].atype = 'C1 ') and (bond^[get_bond(a_view,a_ref)].btype = 'T') and
     (atom^[a_ref].atype = 'N1 ') and (atom^[a_ref].neighbor_count = 2) and
     (atom^[a_view].neighbor_count = 1) then r := true;
  is_isonitrile := r;
end;


function is_cyanate(a_view,a_ref:integer):boolean;
// a_view = C, a_ref = N
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  if is_cyano(a_view,a_ref) then
  begin
    if (atom^[a_view].neighbor_count = 2) then
      begin
        nb := get_nextneighbors(a_view,a_ref);
        if (is_alkoxy(a_view,nb[1])) or
           (is_aryloxy(a_view,nb[1])) then r := true;
      end;
    end;
  is_cyanate := r;
end;


function is_thiocyanate(a_view,a_ref:integer):boolean;
var
  r  : boolean;
  nb : neighbor_rec;
begin
  r := false;
  if is_cyano(a_view,a_ref) then
  begin
    if (atom^[a_view].neighbor_count = 2) then 
      begin
        nb := get_nextneighbors(a_view,a_ref);
        if (is_alkylsulfanyl(a_view,nb[1])) or 
           (is_arylsulfanyl(a_view,nb[1])) then r := true;
      end;
    end;
  is_thiocyanate := r;
end;


procedure update_Htotal;
var
  i, j, b_id : integer;
  nb : neighbor_rec;
  single_count, double_count, triple_count, arom_count : integer;
  total_bonds : integer;
  Htotal : integer;
  nval   : integer;   // new in v0.3
  diazon : boolean;       // new in v0.3j
  nb2    : neighbor_rec;  // new in v0.3j
  a1, a2, a3 : integer;   // new in v0.3j
begin
  if n_atoms < 1 then exit;
  diazon := false;
  fillchar(nb,sizeof(neighbor_rec),0);
  for i := 1 to n_atoms do
    begin
      single_count := 0;
      double_count := 0;
      triple_count := 0;
      arom_count   := 0;
      total_bonds  := 0;
      Htotal    := 0;
      nb := get_neighbors(i);
      if atom^[i].neighbor_count > 0 then
        begin  // count single, double, triple, and aromatic bonds to all neighbor atoms
          for j := 1 to atom^[i].neighbor_count do
            begin
              b_id := get_bond(i,nb[j]);
              if b_id > 0 then
                begin
                  if bond^[b_id].btype = 'S' then inc(single_count);
                  if bond^[b_id].btype = 'D' then inc(double_count);
                  if bond^[b_id].btype = 'T' then inc(triple_count);
                  if bond^[b_id].btype = 'A' then inc(arom_count);
                end;
            end;
          //check for diazonium salts
          a1 := i; a2 := nb[1];
          if (atom^[a1].element = 'N ') and (atom^[a2].element = 'N ') then
            begin
              if (atom^[a2].neighbor_count = 2) then
                begin
                  nb2 := get_nextneighbors(a2,a1);
                  a3 := nb2[1];
                  if (atom^[a3].element = 'C ') and is_diazonium(a3,a2) then diazon := true;
                end;
            end;
        end;
      total_bonds := single_count + 2*double_count + 3*triple_count + trunc(1.5*arom_count);  
      // calculate number of total hydrogens per atom
      //nval := nvalences(atom^[i].element);    // new in v0.3
      nval := atom^[i].nvalences;    // new in v0.3m
      if (atom^[i].element = 'P ') then 
        begin
          if ((total_bonds - atom^[i].formal_charge) > 3) then nval := 5;  // refined in v0.3n
        end;                                  // 
      if (atom^[i].element = 'S ') then       // v0.3h
        begin
          if (total_bonds > 2) and (atom^[i].formal_charge < 1) then nval := 4;  // updated in v0.3j
          if total_bonds > 4 then nval := 6;  // this will need some refinement...
        end;                                  // 
      Htotal := nval - total_bonds + atom^[i].formal_charge;
      if (atom^[i].radical_type = 1) or (atom^[i].radical_type = 3) then Htotal := Htotal - 2; // v0.3p
      if (atom^[i].radical_type = 2) then Htotal := Htotal - 1; // v0.3p
      if diazon then Htotal := 0;      // v0.3j
      if Htotal < 0 then Htotal := 0;  // e.g., N in nitro group
      atom^[i].Htot := Htotal;
      if atom^[i].Hexp > atom^[i].Htot then atom^[i].Htot := atom^[i].Hexp;  // v0.3n; just to be sure...
    end;
end;


procedure update_atypes;
var
  i, j, b_id : integer;
  nb : neighbor_rec;
  single_count, double_count, triple_count, arom_count, acyl_count : integer;
  C_count, O_count : integer;
  total_bonds : integer;
  NdO_count : integer;
  NdC_count : integer;
  Htotal : integer;
begin
  if n_atoms < 1 then exit;
  fillchar(nb,sizeof(neighbor_rec),0);
  for i := 1 to n_atoms do
    begin
      single_count := 0;
      double_count := 0;
      triple_count := 0;
      arom_count   := 0;
      total_bonds  := 0;
      acyl_count   := 0;
      C_count      := 0;
      O_count      := 0;
      NdO_count := 0;
      NdC_count := 0;
      Htotal    := 0;
      nb := get_neighbors(i);
      if atom^[i].neighbor_count > 0 then
        begin  // count single, double, triple, and aromatic bonds to all neighbor atoms
          for j := 1 to atom^[i].neighbor_count do
            begin
              if (is_oxo_C(nb[j])) or (is_thioxo_C(nb[j])) then inc(acyl_count);
              if atom^[(nb[j])].element = 'C ' then inc(C_count);
              if atom^[(nb[j])].element = 'O ' then inc(O_count);
              b_id := get_bond(i,nb[j]);
              if b_id > 0 then
                begin
                  if bond^[b_id].btype = 'S' then inc(single_count);
                  if bond^[b_id].btype = 'D' then inc(double_count);
                  if bond^[b_id].btype = 'T' then inc(triple_count);
                  if bond^[b_id].btype = 'A' then     // v0.3n: special treatment for acyclic bonds
                    begin                             // flagged as "aromatic" (in query structures)
                      if (bond^[b_id].ring_count > 0) then inc(arom_count) else inc(double_count);
                    end;
                  if ((atom^[i].element = 'N ') and (atom^[(nb[j])].element = 'O ')) or
                     ((atom^[i].element = 'O ') and (atom^[(nb[j])].element = 'N ')) then
                     begin
                       // check if it is an N-oxide drawn with a double bond ==> should be N3
                       if bond^[b_id].btype = 'D' then inc(NdO_count);
                     end;
                  if ((atom^[i].element = 'N ') and (atom^[(nb[j])].element = 'C ')) or
                     ((atom^[i].element = 'C ') and (atom^[(nb[j])].element = 'N ')) then
                     begin
                       if bond^[b_id].btype = 'D' then inc(NdC_count);
                     end;
                end;
            end;
          total_bonds := single_count + 2*double_count + 3*triple_count + trunc(1.5*arom_count);  
          // calculate number of total hydrogens per atom
          //Htotal := nvalences(atom^[i].element) - total_bonds + atom^[i].formal_charge;
          Htotal := atom^[i].nvalences - total_bonds + atom^[i].formal_charge;
          if Htotal < 0 then Htotal := 0;  // e.g., N in nitro group
          atom^[i].Htot := Htotal;
          // refine atom types, based on bond types
          if atom^[i].element = 'C ' then
            begin
              if (arom_count > 1) then atom^[i].atype := 'CAR';
              if (triple_count = 1) or (double_count = 2) then atom^[i].atype := 'C1 ';
              if (double_count = 1) then atom^[i].atype := 'C2 ';
              if (triple_count = 0) and (double_count = 0) and (arom_count < 2) then atom^[i].atype := 'C3 ';
            end;  
          if atom^[i].element = 'O ' then
            begin
              if (double_count = 1) then atom^[i].atype := 'O2 ';
              if (double_count = 0) then atom^[i].atype := 'O3 ';
            end;
          if atom^[i].element = 'N ' then
            begin
              if total_bonds > 3 then
                begin
                  if O_count = 0 then
                    begin
                      if (single_count > 3) or
                        ((single_count = 2) and (double_count = 1) and (C_count >=2)) then
                        atom^[i].formal_charge := 1;
                    end else  // could be an N-oxide -> should be found elsewhere 
                    begin
                      if (O_count = 1) and (atom^[i].formal_charge = 0) then 
                        atom^[i].atype := 'N3 ';  // v0.3m
                      if (O_count = 2) and (atom^[i].formal_charge = 0) then 
                        begin
                          if (atom^[i].neighbor_count > 2) then atom^[i].atype := 'N2 ';  // nitro v0.3o
                          if (atom^[i].neighbor_count = 2) then atom^[i].atype := 'N1 ';  // NO2   v0.3o
                        end;
                      // the rest is left empty, so far....
                    end;
                end;
              if (triple_count = 1) or ((double_count = 2) and
                 (atom^[i].neighbor_count = 2)) then atom^[i].atype := 'N1 '; // v0.3n
              if (double_count = 1) then 
                begin
                  //if NdC_count > 0 then atom^[i].atype := 'N2 ';
                  if (NdC_count = 0) and (NdO_count > 0) and
                     (C_count >= 2) then atom^[i].atype := 'N3 '  // N-oxide is N3 except in hetarene etc.
                  else atom^[i].atype := 'N2 ';                   // fallback, added in v0.3g 
                end;  
              if (arom_count > 1) or (atom^[i].arom = true) then atom^[i].atype := 'NAR';  // v0.3n
              if (triple_count = 0) and (double_count = 0) then 
                begin
                  if (atom^[i].formal_charge = 0) then 
                    begin
                      if (acyl_count = 0) then atom^[i].atype := 'N3 ';
                      if (acyl_count > 0) then atom^[i].atype := 'NAM';
                    end;  
                  if (atom^[i].formal_charge = 1) then atom^[i].atype := 'N3+';
                end;
            end;  
          if atom^[i].element = 'P ' then
            begin
              if (single_count > 4) then atom^[i].atype := 'P4 ';
              if (single_count <= 4) and (double_count = 0) then atom^[i].atype := 'P3 ';
              if (double_count = 2) then atom^[i].atype := 'P3D';
            end;
          if atom^[i].element = 'S ' then
            begin
              if (double_count = 1) and (single_count = 0) then atom^[i].atype := 'S2 ';
              if (double_count = 0) then atom^[i].atype := 'S3 ';
              if (double_count = 1) and (single_count > 0) then atom^[i].atype := 'SO ';
              if (double_count = 2) and (single_count > 0) then atom^[i].atype := 'SO2';
            end;
          // further atom types should go here
        end;
    end;
end;


procedure update_atypes_quick;  // v0.4b
var
  i : integer;
begin
  if n_atoms < 1 then exit;
  for i := 1 to n_atoms do
    begin
      if (atom^[i].element = 'C ') and atom^[i].arom then atom^[i].atype := 'CAR';
      if (atom^[i].atype = 'N2 ') and atom^[i].arom then atom^[i].atype := 'NAR';
      if atom^[i].Hexp > atom^[i].Htot then atom^[i].Htot := atom^[i].Hexp;  
    end;
end;


procedure chk_arom;
var
  i, j, pi_count, ring_size : integer;
  b, a1, a2 : integer;  // v0.3n
  testring : ringpath_type;
  a_ref, a_prev, a_next : integer;
  b_bk, b_fw, b_exo : integer;
  bt_bk, bt_fw : char;
  ar_bk, ar_fw, ar_exo : boolean;  // new in v0.3
  conj_intr, ko, aromatic : boolean;
  aromatic_bt : boolean;  // v0.3n
  n_db, n_sb, n_ar : integer;
  cumul : boolean;
  exo_mC : integer;        // v0.3j
  arom_pi_diff : integer;  // v0.3j
begin
  if n_rings < 1 then exit;
  // first, do a very quick check for benzene, pyridine, etc.
  for i := 1 to n_rings do
    begin
      ring_size := ringprop^[i].size;
      if (ring_size = 6) then
        begin
          fillchar(testring,sizeof(ringpath_type),0);
          for j := 1 to ring_size do testring[j] := ring^[i,j];
          cumul := false;
          n_sb := 0;
          n_db := 0;
          n_ar := 0;
          a_prev := testring[ring_size];
          for j := 1 to ring_size do
            begin
              a_ref := testring[j];
              if (j < ring_size) then a_next := testring[(j+1)] else a_next := testring[1];
              b_bk  := get_bond(a_prev,a_ref);
              b_fw  := get_bond(a_ref,a_next);
              bt_bk := bond^[b_bk].btype;
              bt_fw := bond^[b_fw].btype;
              if (bt_fw = 'S') then inc(n_sb);
              if (bt_fw = 'D') then inc(n_db);
              if (bt_fw = 'A') then inc(n_ar);
              if (bt_fw <> 'A') and (bt_bk = bt_fw) then cumul := true;
              a_prev := a_ref;
            end;
          if (n_ar = 6) or ((n_sb = 3) and (n_db = 3) and (cumul = false)) then
            begin   // this ring is aromatic
              a_prev := testring[ring_size];
              for j := 1 to ring_size do
                begin
                  a_ref := testring[j];
                  b_bk  := get_bond(a_prev,a_ref);
                  bond^[b_bk].arom := true;
                  a_prev := a_ref;
                end;
              ringprop^[i].arom := true;            
            end;
        end;
    end;  
  for i := 1 to n_rings do
    begin
      if (ringprop^[i].arom = false) then   
        begin   // do the hard work only for those rings which are not yet flagged aromatic
          fillchar(testring,sizeof(ringpath_type),0);
          ring_size := ringprop^[i].size;  // v0.3j
          for j := 1 to ring_size do testring[j] := ring^[i,j];  // v0.3j
          pi_count  := 0;
          arom_pi_diff := 0;  // v0.3j
          conj_intr := false;
          ko        := false;
          a_prev    := testring[ring_size];
          for j := 1 to ring_size do
            begin
              a_ref := testring[j];
              if (j < ring_size) then a_next := testring[(j+1)] else a_next := testring[1];
              b_bk  := get_bond(a_prev,a_ref);
              b_fw  := get_bond(a_ref,a_next);
              bt_bk := bond^[b_bk].btype;
              bt_fw := bond^[b_fw].btype;
              ar_bk := bond^[b_bk].arom;
              ar_fw := bond^[b_fw].arom;
              if ((bt_bk = 'S') and (bt_fw = 'S') and (ar_bk = false) and (ar_fw = false)) then
                begin
                  // first, assume the worst case (interrupted conjugation)
                  conj_intr := true;  
                  // conjugation can be restored by hetero atoms
                  if (atom^[a_ref].atype = 'O3 ') or (atom^[a_ref].atype = 'S3 ') or
                     (atom^[a_ref].element = 'N ') or (atom^[a_ref].element = 'SE') or
                     (atom^[a_ref].element = 'Q ') then    // v0.4: query atom of type Q is OK
                     begin
                       conj_intr := false;
                       inc(pi_count,2);  // lone pair adds for 2 pi electrons
                     end;
                  // conjugation can be restored by a formal charge at a methylene group
                  if (atom^[a_ref].element = 'C ') and (atom^[a_ref].formal_charge <> 0) then
                    begin
                      conj_intr := false;
                      pi_count  := pi_count - atom^[a_ref].formal_charge;  // neg. charge increases pi_count!
                    end;
                  // conjugation can be restored by carbonyl groups etc.
                  if (is_oxo_C(a_ref)) or (is_thioxo_C(a_ref)) or (is_exocyclic_imino_C(a_ref,i)) then
                    begin
                      conj_intr := false;
                    end;
                  // conjugation can be restored by exocyclic C=C double bond,
                  // adds 2 pi electrons to 5-membered rings, not to 7-membered rings (CAUTION!)
                  // apply only to non-aromatic exocyclic C=C bonds
                  exo_mC := find_exocyclic_methylene_C(a_ref,i);  // v0.3j
                  if ((exo_mC > 0) and odd(ring_size)) then       // v0.3j
                    begin
                      b_exo  := get_bond(a_ref,exo_mC);           // v0.3j 
                      ar_exo := bond^[b_exo].arom;
                      if ((ring_size - 1) mod 4 = 0) then  // 5-membered rings and related
                        begin
                          conj_intr := false;
                          inc(pi_count,2);
                        end else                           // 7-membered rings and related
                        begin
                          if not ar_exo then conj_intr := false;
                        end;
                    end;
                  // if conjugation is still interrupted ==> knock-out
                  if conj_intr then ko := true;
                end else
                begin
                  if ((bt_bk = 'S') and (bt_fw = 'S') and (ar_bk = true) and (ar_fw = true)) then
                    begin
                      if (atom^[a_ref].atype = 'O3 ') or (atom^[a_ref].atype = 'S3 ') or
                         (atom^[a_ref].element = 'N ') or (atom^[a_ref].element = 'SE') then
                         begin
                           inc(pi_count,2);  // lone pair adds for 2 pi electrons
                         end;
                      if (atom^[a_ref].element = 'C ') and (atom^[a_ref].formal_charge <> 0) then
                        begin
                          pi_count  := pi_count - atom^[a_ref].formal_charge;  // neg. charge increases pi_count!
                        end;
                      exo_mC := find_exocyclic_methylene_C(a_ref,i);  // v0.3j
                      if ((exo_mC > 0) and odd(ring_size)) then       // v0.3j
                        begin
                          b_exo := get_bond(a_ref,exo_mC);            // v0.3j
                          ar_exo := bond^[b_exo].arom;
                          if ((ring_size - 1) mod 4 = 0) then  // 5-membered rings and related
                            begin
                              inc(pi_count,2);
                            end;
                        end;
                    end else    // any other case: increase pi count by one electron
                    begin
                      inc(pi_count);  // v0.3j; adjustment for bridgehead N: see below
                      if ((bt_bk = 'S') and (bt_fw = 'S')) and
                          (((ar_bk = true) and (ar_fw = false)) or
                           ((ar_bk = false) and (ar_fw = true)) ) then
                        begin
                          // v0.3j; if a bridgehead N were not aromatic, it could 
                          // contribute 2 pi electrons --> try also this variant
                          // (example: CAS 32278-54-9)
                          if (atom^[a_ref].element = 'N ') then inc(arom_pi_diff);
                        end;
                    end;
                end;
              // last command:
              a_prev := a_ref;
            end;  // for j := 1 to ring_size
          // now we can draw our conclusion
          //if not ((ko) or (odd(pi_count))) then
          if not ko then    // v0.3j; odd pi_count might be compensated by arom_pi_diff
            begin  // apply Hueckel's rule
              if (abs(ring_size - pi_count) < 2) and 
                 (((pi_count - 2) mod 4 = 0) or
                  (((pi_count + arom_pi_diff) - 2) mod 4 = 0) ) then
                begin
                  // this ring is aromatic
                  ringprop^[i].arom := true;
                  // now mark _all_ bonds in the ring as aromatic
                  a_prev := testring[ring_size];
                  for j := 1 to ring_size do
                    begin
                      a_ref := testring[j];
                      bond^[get_bond(a_prev,a_ref)].arom := true;
                      a_prev := a_ref;
                     end;
                end;
            end;
        end;
    end;  // (for i := 1 to n_rings)
  // finally, mark all involved atoms as aromatic
  for i := 1 to n_bonds do
    begin
      if bond^[i].arom then
        begin
          a1 := bond^[i].a1;  // v0.3n
          a2 := bond^[i].a2;  // v0.3n
          atom^[a1].arom := true;
          atom^[a2].arom := true;
          // v0.3n: update atom types if applicable (C and N)
          if (atom^[a1].element = 'C ') then atom^[a1].atype := 'CAR';
          if (atom^[a2].element = 'C ') then atom^[a2].atype := 'CAR';
          if (atom^[a1].element = 'N ') then atom^[a1].atype := 'NAR';
          if (atom^[a2].element = 'N ') then atom^[a2].atype := 'NAR';
        end;
    end;
  // update aromaticity information in ringprop
  // new in v0.3n: accept rings as aromatic if all bonds are of type 'A'
  for i := 1 to n_rings do
    begin
      testring := ring^[i];
      //ring_size := path_length(testring);
      ring_size := ringprop^[i].size;  // v0.3j
      aromatic := true;
      aromatic_bt := true;  // v0.3n
      a_prev := testring[ring_size];
      for j := 1 to ring_size do
        begin
          a_ref := testring[j];
          b := get_bond(a_prev,a_ref);  // v0.3n
          if (not (bond^[b].arom)) then aromatic := false;           
          if (not (bond^[b].btype = 'A')) then aromatic_bt := false;  // v0.3n
          a_prev := a_ref;
        end;
      if aromatic_bt and (not aromatic) then  // v0.3n: update aromaticity flag
        begin
          a_prev := testring[ring_size];
          for j := 1 to ring_size do
            begin
              a_ref := testring[j];
              b := get_bond(a_prev,a_ref);
              bond^[b].arom := true;
              if (atom^[a_ref].element = 'C ') then atom^[a_ref].atype := 'CAR';
              if (atom^[a_ref].element = 'N ') then atom^[a_ref].atype := 'NAR';
              a_prev := a_ref;
            end;
          aromatic := true;
        end;                                  // end v0.3n block  
      if aromatic then ringprop^[i].arom := true else ringprop^[i].arom := false;
    end;  
end;


procedure write_mol;
var
  i, j : integer;
  testring : ringpath_type;
  ring_size : integer;
  //aromatic : boolean;
  //a_prev, a_ref : integer;
begin
  if progmode = pmCheckMol then
    writeln('Molecule name: ',molname)
  else
    writeln('Molecule name (haystack): ',molname);
  writeln('atoms: ',n_atoms,'  bonds: ',n_bonds,'  rings: ',n_rings);
  if n_atoms < 1 then exit;
  if n_bonds < 1 then exit;
  for i := 1 to n_atoms do
    begin
      if i <   10 then write(' ');
      if i <  100 then write(' ');
      if i < 1000 then write(' ');
      write(i,' ',atom^[i].element,' ',atom^[i].atype,' ',atom^[i].x:9:4,' ',atom^[i].y:9:4,' ');
      write(atom^[i].z:9:4);
      write('  (',atom^[i].neighbor_count,' heavy-atom neighbors, Hexp: ',atom^[i].Hexp,' Htot: ',atom^[i].Htot,')');
      if (atom^[i].formal_charge <> 0) then write('  charge: ',atom^[i].formal_charge);
      if (atom^[i].arom) then write(' aromatic');  // 
      writeln;
    end;
  for i := 1 to n_bonds do
    begin
      if i <   10 then write(' ');
      if i <  100 then write(' ');
      if i < 1000 then write(' ');
      write(i,' ',bond^[i].a1,' ',bond^[i].a2,' ',bond^[i].btype);
      if bond^[i].ring_count > 0 then write(', contained in ',bond^[i].ring_count,' ring(s)');
      if bond^[i].arom then write(' (aromatic) ');
      writeln;
    end;
  if n_rings > 0 then
    begin
      for i := 1 to n_rings do
        begin
          write('ring ',i,': ');
          //aromatic := true;
          fillchar(testring,sizeof(ringpath_type),0);
          ring_size := ringprop^[i].size;  // v0.3j
          //for j := 1 to max_ringsize do if ring^[i,j] > 0 then testring[j] := ring^[i,j];
          for j := 1 to ring_size do testring[j] := ring^[i,j];  // v0.3j
          //ring_size := path_length(testring);
          //a_prev := testring[ring_size];
          for j := 1 to ring_size do
            begin
              write(testring[j],' ');
              //a_ref := testring[j];
              //if (not bond^[get_bond(a_prev,a_ref)].arom) then aromatic := false;
              //a_prev := a_ref;
            end;
          //if aromatic then write(' (aromatic)');
          if (ringprop^[i].arom) then write(' (aromatic)');
          if (ringprop^[i].envelope) then write(' (env)');
          writeln;
        end;
    end;
end;


procedure write_needle_mol;
var
  i, j : integer;
  testring : ringpath_type;
  ring_size : integer;
  aromatic : boolean;
  a_prev, a_ref : integer;
begin
  writeln('Molecule name (needle): ',ndl_molname);
  writeln('atoms: ',ndl_n_atoms,'  bonds: ',ndl_n_bonds,'  rings: ',ndl_n_rings);
  if ndl_n_atoms < 1 then exit;
  if ndl_n_bonds < 1 then exit;
  for i := 1 to ndl_n_atoms do
    begin
      if i <   10 then write(' ');
      if i <  100 then write(' ');
      if i < 1000 then write(' ');
      write(i,' ',ndl_atom^[i].element,' ',ndl_atom^[i].atype,' ',ndl_atom^[i].x:9:4,' ',atom^[i].y:9:4,' ');
      write(ndl_atom^[i].z:9:4);
      write('  (',ndl_atom^[i].neighbor_count,' heavy-atom neighbors, Hexp: ',ndl_atom^[i].Hexp,' Htot: ',ndl_atom^[i].Htot,')');
      if (ndl_atom^[i].formal_charge <> 0) then write('  charge: ',ndl_atom^[i].formal_charge);
      writeln;
    end;
  for i := 1 to ndl_n_bonds do
    begin
      if i <   10 then write(' ');
      if i <  100 then write(' ');
      if i < 1000 then write(' ');
      write(i,' ',ndl_bond^[i].a1,' ',ndl_bond^[i].a2,' ',ndl_bond^[i].btype);
      if ndl_bond^[i].ring_count > 0 then write(', contained in ',ndl_bond^[i].ring_count,' ring(s)');
      if ndl_bond^[i].arom then write(' (aromatic) ');
      writeln;
    end;
  if ndl_n_rings > 0 then
    begin
      for i := 1 to ndl_n_rings do
        begin
          aromatic := true;
          fillchar(testring,sizeof(ringpath_type),0);
          for j := 1 to max_ringsize do if ndl_ring^[i,j] > 0 then testring[j] := ndl_ring^[i,j];
          ring_size := path_length(testring);
          write('ring ',i,': ');
          a_prev := testring[ring_size];
          for j := 1 to ring_size do
            begin
              write(testring[j],' ');
              a_ref := testring[j];
              if (not ndl_bond^[get_ndl_bond(a_prev,a_ref)].arom) then aromatic := false;  // v0.3k
              a_prev := a_ref;
            end;
          if aromatic then write(' (aromatic)');
          writeln;
        end;
    end;
end;


procedure chk_so2_deriv(a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  het_count, o_count, or_count, hal_count, n_count, c_count : integer;
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  if atom^[a_ref].atype <>'SO2' then exit;
  nb := get_neighbors(a_ref);
  het_count := 0; o_count := 0; or_count := 0; hal_count := 0; n_count := 0; c_count := 0;
  for i := 1 to atom^[a_ref].neighbor_count do
    begin
      if bond^[(get_bond(a_ref,nb[i]))].btype = 'S' then
        begin
          nb_el := atom^[(nb[i])].element;
          if (nb_el <> 'C ') and (nb_el <> 'H ') and (nb_el <> 'D ') and  // added 'D ' in v0.3n
             (nb_el <> 'DU') and (nb_el <> 'LP') then inc(het_count);
          if nb_el = 'O ' then
            begin
              inc(o_count);
              if is_alkoxy(a_ref,nb[i]) or is_aryloxy(a_ref,nb[i]) then inc(or_count);
            end;
          if nb_el = 'N ' then inc(n_count);
          if nb_el = 'C ' then inc(c_count);
          if (nb_el = 'F ') or (nb_el = 'CL') or (nb_el = 'BR') or (nb_el = 'I ') then inc(hal_count);
        end;
    end;
  if het_count = 2 then   // sulfuric acid derivative
    begin
      fg[fg_sulfuric_acid_deriv] := true;
      if opt_pos then add2fgloc(fg_sulfuric_acid_deriv,a_ref);   // v0.5
      if (o_count = 2) then
        begin
          if (or_count = 0) then 
            begin
              fg[fg_sulfuric_acid] := true;
              if opt_pos then add2fgloc(fg_sulfuric_acid,a_ref);   // v0.5
            end;
          if (or_count = 1) then 
            begin
              fg[fg_sulfuric_acid_monoester] := true;
              if opt_pos then add2fgloc(fg_sulfuric_acid_monoester,a_ref);   // v0.5
            end;
          if (or_count = 2) then 
            begin
              fg[fg_sulfuric_acid_diester] := true;
              if opt_pos then add2fgloc(fg_sulfuric_acid_diester,a_ref);   // v0.5
            end;
        end;
      if (o_count = 1) then
        begin
          if (or_count = 1) and (n_count = 1) then 
            begin
              fg[fg_sulfuric_acid_amide_ester] := true;
              if opt_pos then add2fgloc(fg_sulfuric_acid_amide_ester,a_ref);   // v0.5
            end;
          if (or_count = 0) and (n_count = 1) then 
            begin
              fg[fg_sulfuric_acid_amide] := true;
              if opt_pos then add2fgloc(fg_sulfuric_acid_amide,a_ref);   // v0.5
            end;
        end;
      if (n_count = 2)   then 
        begin
          fg[fg_sulfuric_acid_diamide] := true;
          if opt_pos then add2fgloc(fg_sulfuric_acid_diamide,a_ref);   // v0.5
        end;
      if (hal_count > 0) then 
        begin
          fg[fg_sulfuryl_halide] := true;
          if opt_pos then add2fgloc(fg_sulfuryl_halide,a_ref);   // v0.5
        end;
    end;
  if (het_count = 1) and (c_count = 1) then   // sulfonic acid derivative
    begin
      fg[fg_sulfonic_acid_deriv] := true;
      if opt_pos then add2fgloc(fg_sulfonic_acid_deriv,a_ref);   // v0.5
      if (o_count = 1) and (or_count = 0) then 
        begin
          fg[fg_sulfonic_acid] := true;
          if opt_pos then add2fgloc(fg_sulfonic_acid,a_ref);   // v0.5
        end;
      if (o_count = 1) and (or_count = 1) then 
        begin
          fg[fg_sulfonic_acid_ester] := true;
          if opt_pos then add2fgloc(fg_sulfonic_acid_ester,a_ref);   // v0.5
        end;
      if (n_count = 1) then 
        begin
          fg[fg_sulfonamide] := true;
          if opt_pos then add2fgloc(fg_sulfonamide,a_ref);   // v0.5
        end;
      if (hal_count = 1) then 
        begin
          fg[fg_sulfonyl_halide] := true;
          if opt_pos then add2fgloc(fg_sulfonyl_halide,a_ref);   // v0.5
        end;
    end;
  if (het_count = 0) and (c_count = 2) then   // sulfone
    begin
      fg[fg_sulfone] := true;
      if opt_pos then add2fgloc(fg_sulfone,a_ref);   // v0.5
    end;
end;


procedure chk_p_deriv(a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  nb_el, dbl_het : str2;
  het_count, oh_count, or_count, hal_count, n_count, c_count : integer;
begin
  if atom^[a_ref].element <>'P ' then exit;
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_ref);
  dbl_het := '';
  het_count := 0; oh_count := 0; or_count := 0; hal_count := 0; n_count := 0; c_count := 0;
  for i := 1 to atom^[a_ref].neighbor_count do
    begin
      if bond^[(get_bond(a_ref,nb[i]))].btype = 'D' then
        begin
          dbl_het := atom^[(nb[i])].element;
        end;
      if bond^[(get_bond(a_ref,nb[i]))].btype = 'S' then
        begin
          nb_el := atom^[(nb[i])].element;
          if (nb_el = 'C ') then inc(c_count);
          if (is_hydroxy(a_ref,nb[i])) then inc(oh_count);
          if (is_alkoxy(a_ref,nb[i])) or (is_aryloxy(a_ref,nb[i])) then inc(or_count);
          if (nb_el = 'N ') then inc(n_count);
          if (nb_el = 'F ') or (nb_el = 'CL') or (nb_el = 'BR') or (nb_el = 'I ') then inc(hal_count);
        end;
    end;
  het_count := oh_count + or_count + hal_count + n_count;
  if (atom^[a_ref].atype = 'P3D') or (atom^[a_ref].atype = 'P4 ') then
    begin
      if (dbl_het = 'O ') then
        begin
          if (c_count = 0) then
            begin
              fg[fg_phosphoric_acid_deriv] := true;
              if opt_pos then add2fgloc(fg_phosphoric_acid_deriv,a_ref);  // v0.5
              if (oh_count = 3)   then 
                begin
                  fg[fg_phosphoric_acid]        := true;
                  if opt_pos then add2fgloc(fg_phosphoric_acid,a_ref);  // v0.5
                end;
              if (or_count > 0)   then 
                begin
                  fg[fg_phosphoric_acid_ester]  := true;
                  if opt_pos then add2fgloc(fg_phosphoric_acid_ester,a_ref);  // v0.5
                end;
              if (hal_count > 0)  then 
                begin
                  fg[fg_phosphoric_acid_halide] := true;            
                  if opt_pos then add2fgloc(fg_phosphoric_acid_halide,a_ref);  // v0.5
                end;
              if (n_count > 0)    then 
                begin
                  fg[fg_phosphoric_acid_amide]  := true;
                  if opt_pos then add2fgloc(fg_phosphoric_acid_amide,a_ref);  // v0.5
                end;
            end;
          if (c_count = 1) then
            begin
              fg[fg_phosphonic_acid_deriv] := true;
              if opt_pos then add2fgloc(fg_phosphonic_acid_deriv,a_ref);  // v0.5
              if (oh_count = 2)   then 
                begin
                  fg[fg_phosphonic_acid]        := true;
                  if opt_pos then add2fgloc(fg_phosphonic_acid,a_ref);  // v0.5
                end;
              if (or_count > 0)   then 
                begin
                  fg[fg_phosphonic_acid_ester]  := true;
                  if opt_pos then add2fgloc(fg_phosphonic_acid_ester,a_ref);  // v0.5
                end;
              //if (hal_count > 0)  then fg[fg_phosphonic_acid_halide] := true;            
              //if (n_count > 0)    then fg[fg_phosphonic_acid_amide]  := true;
            end;
          if (c_count = 3) then 
            begin
              fg[fg_phosphinoxide] := true;  
              if opt_pos then add2fgloc(fg_phosphinoxide,a_ref);  // v0.5
            end;
        end;
      if (dbl_het = 'S ') then
        begin
          if (c_count = 0) then
            begin
              fg[fg_thiophosphoric_acid_deriv] := true;
              if opt_pos then add2fgloc(fg_thiophosphoric_acid_deriv,a_ref);  // v0.5
              if (oh_count = 3)   then 
                begin
                  fg[fg_thiophosphoric_acid]        := true;
                  if opt_pos then add2fgloc(fg_thiophosphoric_acid,a_ref);  // v0.5
                end;
              if (or_count > 0)   then 
                begin
                  fg[fg_thiophosphoric_acid_ester]  := true;
                  if opt_pos then add2fgloc(fg_thiophosphoric_acid_ester,a_ref);  // v0.5
                end;
              if (hal_count > 0)  then 
                begin
                  fg[fg_thiophosphoric_acid_halide] := true;            
                  if opt_pos then add2fgloc(fg_thiophosphoric_acid_halide,a_ref);  // v0.5
                end;
              if (n_count > 0)    then 
                begin
                  fg[fg_thiophosphoric_acid_amide]  := true;
                  if opt_pos then add2fgloc(fg_thiophosphoric_acid_amide,a_ref);  // v0.5
                end;
            end;
        end;
    end;
//  if (atom^[a_ref].atype = 'P4 ') then fg[fg_phosphoric_acid_deriv] := true;
  if (atom^[a_ref].atype = 'P3 ') then    // changed P3D into P3 in v0.3b
    begin
      if (c_count = 3) and (het_count = 0) then 
        begin
          fg[fg_phosphine] := true;
          if opt_pos then add2fgloc(fg_phosphine,a_ref);  // v0.5
        end;
      if (c_count = 3) and (oh_count = 1) then 
        begin
          fg[fg_phosphinoxide] := true;
          if opt_pos then add2fgloc(fg_phosphinoxide,a_ref);  // v0.5
        end;
    end;
end;


procedure chk_b_deriv(a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  het_count, oh_count, or_count, hal_count, n_count, c_count : integer;
begin
  if atom^[a_ref].element <>'B ' then exit;
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_ref);
  het_count := 0; oh_count := 0; or_count := 0; hal_count := 0; n_count := 0; c_count := 0;
  for i := 1 to atom^[a_ref].neighbor_count do
    begin
      if bond^[(get_bond(a_ref,nb[i]))].btype = 'S' then
        begin
          nb_el := atom^[(nb[i])].element;
          if (nb_el = 'C ') then inc(c_count) else if ((nb_el <> 'H ') and (nb_el <> 'D ') and // v0.3n: D
                                                      (nb_el <>'LP')) then inc(het_count);
          if (is_hydroxy(a_ref,nb[i])) then inc(oh_count);
          if (is_alkoxy(a_ref,nb[i])) or (is_aryloxy(a_ref,nb[i])) then inc(or_count);  // fixed in v0.3b
          if (nb_el = 'N ') then inc(n_count);
          if (nb_el = 'F ') or (nb_el = 'CL') or (nb_el = 'BR') or (nb_el = 'I ') then inc(hal_count);
        end;
    end;
  het_count := oh_count + or_count + hal_count + n_count;  // fixed in v0.3b
  if (c_count = 1) and (het_count = 2) then
    begin
      fg[fg_boronic_acid_deriv] := true;
      if opt_pos then add2fgloc(fg_boronic_acid_deriv,a_ref);  // v0.5
      if (oh_count = 2)   then 
        begin
          fg[fg_boronic_acid]        := true;
          if opt_pos then add2fgloc(fg_boronic_acid,a_ref);  // v0.5
        end;
      if (or_count > 0)   then 
        begin
          fg[fg_boronic_acid_ester]  := true;
          if opt_pos then add2fgloc(fg_boronic_acid_ester,a_ref);  // v0.5
        end;
    end;
end;


procedure chk_ammon(a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  het_count, o_count, or_count, r_count : integer;
  bt : char;  // v0.3k
  bo_sum : single;
  ha : boolean;
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  if (atom^[a_ref].atype <>'N3+') and (atom^[a_ref].formal_charge = 0) then exit;
  if (atom^[a_ref].element <>'N ') then exit;   // just to be sure;  v0.3i
  nb := get_neighbors(a_ref);
  het_count := 0; o_count := 0; or_count := 0; r_count := 0;
  bo_sum := 0;
  for i := 1 to atom^[a_ref].neighbor_count do
    begin
      bt := bond^[(get_bond(a_ref,nb[i]))].btype;  // v0.3k
      nb_el := atom^[(nb[i])].element;             // v0.3k
      ha := atom^[nb[i]].heavy;                   // v0.3k
      if (bt = 'S') then
        begin
          if ha then bo_sum := bo_sum + 1;
          if (nb_el <> 'C ') and (nb_el <> 'H ') and
             (nb_el <> 'D ') and (nb_el <> 'DU') then   // added 'D ' in v0.3n
             begin
               inc(het_count);
               if (nb_el = 'O ') then begin
                 inc(o_count);
                 if atom^[(nb[i])].neighbor_count > 1 then inc(or_count);
               end;  
             end;
          if (is_alkyl(a_ref,nb[i])) or (is_aryl(a_ref,nb[i]))   // v0.3k
              or (is_alkenyl(a_ref,nb[i])) or (is_alkynyl(a_ref,nb[i])) then inc(r_count);
        end; 
      if (bt = 'D') then
        begin
          if ha then bo_sum := bo_sum + 2;
          if (nb_el <> 'C ') then
            begin
              inc(het_count,2);
              if (nb_el = 'O ') then 
                begin
                  inc(o_count,2);
                end;  
            end;
          if nb_el = 'C ' then inc(r_count);
        end; 
      if (bt = 'A') and ha then bo_sum := bo_sum + 1.5;
    end;   // v0.3k: corrected end of "for ..." loop
  if (het_count = 0) and (r_count = 4) then 
    begin
      fg[fg_quart_ammonium] := true;
      if opt_pos then add2fgloc(fg_quart_ammonium,a_ref);   // v0.5
    end;
  if (het_count = 1) and (atom^[a_ref].neighbor_count >= 3) then
    begin
      if ((o_count = 1) and (or_count = 0) and (bo_sum > 3)) then 
        begin
          fg[fg_n_oxide] := true;  // finds only aliphatic N-oxides!
          if opt_pos then add2fgloc(fg_n_oxide,a_ref);   // v0.5
        end;
      if ((((o_count = 1) and (or_count = 1)) or (o_count = 0))
       and (atom^[a_ref].arom = true)) then 
         begin
           fg[fg_quart_ammonium] := true;
           if opt_pos then add2fgloc(fg_quart_ammonium,a_ref);   // v0.5
         end;
    end;
end;


procedure swap_atoms(var a1,a2:integer);
var
  a_tmp : integer;
begin
  a_tmp := a1;
  a1 := a2;
  a2 := a_tmp;
end;


procedure orient_bond(var a1,a2:integer);
var
  a1_el,a2_el : str2;
begin
  a1_el := atom^[a1].element;
  a2_el := atom^[a2].element;
  if (a1_el = 'H ') or (a2_el = 'H ') or  (a1_el = 'D ') or (a2_el = 'D ') then exit; // v0.3n: D
  if (a2_el = 'C ') and (a1_el <> 'C ') then swap_atoms(a1,a2);
  if (a2_el = a1_el) then
    begin
      if hetbond_count(a1) > hetbond_count(a2) then swap_atoms(a1,a2);
    end;
  if (a2_el <> 'C ') and (a1_el <> 'C ') and (a1_el <> a2_el) then
    begin
      if (a1_el = 'O ') or (a2_el = 'O ') then
        begin
          if (a1_el = 'O ') then swap_atoms(a1,a2);
        end;
    end;
  if (a2_el <> 'C ') and (a1_el <> 'C ') and (a1_el = a2_el) then
    begin
      if ((atom^[a2].neighbor_count - hetbond_count(a2)) > 
          (atom^[a1].neighbor_count - hetbond_count(a1))) then swap_atoms(a1,a2);
    end;
end;


procedure chk_imine(a_ref,a_view:integer);
// a_ref = C, a_view = N
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  a_het : integer;
  a_c : integer;
  het_count, c_count : integer;
  o_count : integer;  // v0.3k
begin
  het_count := 0; c_count := 0;
  o_count := 0;  // v0.3k
  if (atom^[a_view].neighbor_count = 1) then
    begin
      if (atom^[a_ref].arom = false) then 
        begin
          fg[fg_imine] := true;
          if opt_pos then add2fgloc(fg_imine,a_ref);  // v0.5
        end;
    end else
    begin
      fillchar(nb,sizeof(neighbor_rec),0);
      nb := get_neighbors(a_view);
      if atom^[a_view].neighbor_count > 1 then
        begin
          for i := 1 to atom^[a_view].neighbor_count do
            begin
              if (nb[i] <> a_ref) and (bond^[(get_bond(a_view,nb[i]))].btype = 'S' ) then
                begin
                  nb_el := atom^[(nb[i])].element;
                  if (nb_el = 'C ') then 
                    begin
                      a_c := nb[i];
                      inc(c_count);
                    end;
                  if (nb_el = 'O ') or (nb_el = 'N ') then
                    begin
                      a_het := nb[i];
                      inc(het_count);
                    end;
                  if (nb_el = 'O ') and (atom^[(nb[i])].neighbor_count = 1) and  // v0.3k
                     (bond^[(get_bond(a_view,nb[i]))].arom = false) then inc(o_count);
                end; 
              if (nb[i] <> a_ref) and (bond^[(get_bond(a_view,nb[i]))].btype = 'D' ) then
                begin    // v0.3k; make sure we do not count nitro groups in "azi" form etc.
                  nb_el := atom^[(nb[i])].element;
                  if (nb_el = 'O ') or (nb_el = 'N ') or (nb_el = 'S ') then 
                    begin
                      a_het := nb[i];  // v0.3m
                      inc(het_count);
                    end;
                  if (nb_el = 'O ') and (atom^[(nb[i])].neighbor_count = 1) and  // v0.3k
                     (bond^[(get_bond(a_view,nb[i]))].arom = false) then inc(o_count);
                end;
            end;
          if (c_count = 1) then 
            begin
              if ((is_alkyl(a_view,a_c)) or (is_aryl(a_view,a_c)) or
                  (is_alkenyl(a_view,a_c)) or (is_alkynyl(a_view,a_c))) and
                 (atom^[a_ref].arom = false) and
                 (het_count = 0) then 
                   begin
                     fg[fg_imine] := true;   // v0.3k
                     if opt_pos then add2fgloc(fg_imine,a_ref);  // v0.5
                   end;
            end;
          if (het_count = 1) then 
            begin
              nb_el := atom^[a_het].element;
              if (nb_el = 'O ') then
                begin
                  if (is_hydroxy(a_view,a_het)) then 
                    begin
                      fg[fg_oxime] := true;
                      if opt_pos then add2fgloc(fg_oxime,a_ref);  // v0.5
                    end;
                  if (is_alkoxy(a_view,a_het))  or (is_aryloxy(a_view,a_het)) or
                     (is_alkenyloxy(a_view,a_het))  or (is_alkynyloxy(a_view,a_het)) then 
                    begin
                      fg[fg_oxime_ether] := true;                  
                      if opt_pos then add2fgloc(fg_oxime_ether,a_ref);  // v0.5
                    end;
                end;
              if (nb_el = 'N ') then
                begin
                  if (is_amino(a_view,a_het)) or (is_alkylamino(a_view,a_het)) or
                     (is_dialkylamino(a_view,a_het)) or (is_alkylarylamino(a_view,a_het)) or
                     (is_arylamino(a_view,a_het)) or (is_diarylamino(a_view,a_het)) then 
                    begin
                      fg[fg_hydrazone] := true;
                      if opt_pos then add2fgloc(fg_hydrazone,a_ref);  // v0.5
                    end else
                    begin  // check for semicarbazone or thiosemicarbazone
                      fillchar(nb,sizeof(neighbor_rec),0);
                      nb := get_neighbors(a_het);
                      if atom^[a_het].neighbor_count > 1 then
                        begin
                          for i := 1 to atom^[a_het].neighbor_count do
                            begin
                              if (nb[i] <> a_view) then
                                begin
                                  if (is_carbamoyl(a_het,nb[i])) then 
                                    begin
                                      fg[fg_semicarbazone] := true;
                                      if opt_pos then add2fgloc(fg_semicarbazone,a_ref);  // v0.5
                                    end;
                                  if (is_thiocarbamoyl(a_het,nb[i])) then 
                                    begin
                                      fg[fg_thiosemicarbazone] := true;                                
                                      if opt_pos then add2fgloc(fg_thiosemicarbazone,a_ref);  // v0.5
                                    end;
                                end;
                            end;
                        end; 
                    end;
                end;
            end;     // v0.3k: nitro groups in "azi" form
          if (het_count = 2) and (o_count = 2) then 
            begin
              fg[fg_nitro_compound] := true;
              if opt_pos then add2fgloc(fg_nitro_compound,a_view);  // v0.5
            end;
        end;
    end;
end;


procedure chk_carbonyl_deriv(a_view,a_ref:integer);
// a_view = C
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  c_count : integer;
  cn_count : integer;
  bt   : char;     // new in v0.3b
  n_db : integer;  // new in v0.3b
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_view);
  c_count := 0; cn_count := 0;
  n_db := 0;  // new in v0.3b
  for i := 1 to atom^[a_view].neighbor_count do
    begin
      bt := bond^[(get_bond(a_view,nb[i]))].btype;
      if bt = 'S' then
        begin
          nb_el := atom^[(nb[i])].element;
          if (nb_el = 'C ') then 
            begin
              if (is_cyano_c(nb[i])) then 
                begin
                  inc(cn_count) 
                end else 
                  begin 
                    inc(c_count);
                  end;
            end;
        end else   // new in v0.3b
        begin
          if bt = 'D' then inc(n_db);
        end;
    end;     
  if is_oxo_c(a_view) then
    begin
      fg[fg_carbonyl]  := true;
      if opt_pos then add2fgloc(fg_carbonyl,a_view);  // v0.5
      if (c_count + cn_count < 2) then 
        begin   // new in v0.3b (detection of ketenes)
          if n_db <= 1 then 
            begin
              fg[fg_aldehyde]  := true;
              if opt_pos then add2fgloc(fg_aldehyde,a_view);  // v0.5
            end else 
            begin
              fg[fg_ketene] := true;
              if opt_pos then add2fgloc(fg_ketene,a_view);  // v0.5
            end;
        end;
      if (c_count = 2) then 
        begin
          if (atom^[a_view].arom) then
            begin
              fg[fg_oxohetarene] := true;
              if opt_pos then add2fgloc(fg_oxohetarene,a_view);  // v0.5
            end else 
            begin
              fg[fg_ketone]    := true;
              if opt_pos then add2fgloc(fg_ketone,a_view);  // v0.5
            end;
        end;
      if (cn_count > 0) then 
        begin
          fg[fg_acyl_cyanide] := true;
          if opt_pos then add2fgloc(fg_acyl_cyanide,a_view);  // v0.5
        end;
    end;
  if is_thioxo_c(a_view) then
    begin
      fg[fg_thiocarbonyl]  := true;
      if opt_pos then add2fgloc(fg_thiocarbonyl,a_view);  // v0.5
      if (c_count < 2) then 
        begin
          fg[fg_thioaldehyde]  := true;
          if opt_pos then add2fgloc(fg_thioaldehyde,a_view);  // v0.5
        end;
      if (c_count = 2) then
        begin
          if (atom^[a_view].arom) then
            begin
              fg[fg_thioxohetarene] := true;
              if opt_pos then add2fgloc(fg_thioxohetarene,a_view);  // v0.5
            end else 
            begin
              fg[fg_thioketone]    := true;
              if opt_pos then add2fgloc(fg_thioketone,a_view);  // v0.5
            end;
        end;
    end;
  if is_imino_c(a_view) then
    begin
      chk_imine(a_view,a_ref);
    end;
end;


procedure chk_carboxyl_deriv(a_view,a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  o_count, n_count, s_count : integer;
  a_o, a_n, a_s : integer;
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_view);
  o_count := 0; n_count := 0; s_count := 0;
  for i := 1 to atom^[a_view].neighbor_count do
    begin
      if bond^[(get_bond(a_view,nb[i]))].btype = 'S' then
        begin
          nb_el := atom^[(nb[i])].element;
          if (nb_el <> 'C ') then 
            begin
              if nb_el = 'O ' then begin inc(o_count); a_o := nb[i]; end;
              if nb_el = 'N ' then begin inc(n_count); a_n := nb[i]; end;
              if nb_el = 'S ' then begin inc(s_count); a_s := nb[i]; end;
            end;
        end;
    end;     
  if (is_oxo_c(a_view)) then
    begin
      if (o_count = 1) then 
        begin  // anhydride is checked somewhere else
          if (bond^[get_bond(a_view,a_o)].arom = false) then 
            begin
              fg[fg_carboxylic_acid_deriv]  := true;
              if opt_pos then add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
            end;
          if (is_hydroxy(a_view,a_o)) then 
            begin
              if (atom^[a_o].formal_charge =  0) then 
                begin
                  fg[fg_carboxylic_acid] := true;
                  if opt_pos then add2fgloc(fg_carboxylic_acid,a_view);   // v0.5
                end;
              if (atom^[a_o].formal_charge = -1) then 
                begin
                  fg[fg_carboxylic_acid_salt] := true;            
                  if opt_pos then add2fgloc(fg_carboxylic_acid_salt,a_view);   // v0.5
                end;
            end;
          if (is_alkoxy(a_view,a_o) or is_aryloxy(a_view,a_o) or
              is_alkenyloxy(a_view,a_o) or is_alkynyloxy(a_view,a_o)) then 
            begin
              if (bond^[get_bond(a_view,a_o)].arom = false) then 
                begin
                  fg[fg_carboxylic_acid_ester] := true;
                  if opt_pos then add2fgloc(fg_carboxylic_acid_ester,a_view);   // v0.5
                end;
              if (bond^[get_bond(a_view,a_o)].ring_count > 0) then
                begin
                  if (bond^[get_bond(a_view,a_o)].arom = true) then
                    begin
                      fg[fg_oxohetarene] := true;
                      if opt_pos then add2fgloc(fg_oxohetarene,a_view);   // v0.5
                    end else 
                    begin
                      fg[fg_lactone] := true;
                      if opt_pos then add2fgloc(fg_lactone,a_view);   // v0.5
                    end;
                end;
            end;            
        end;
      if (n_count = 1) then 
        begin
          if (bond^[get_bond(a_view,a_n)].arom = false) then
            begin 
              fg[fg_carboxylic_acid_deriv]  := true;
              if opt_pos then add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
            end else
            begin
              //fg[fg_lactam_heteroarom] := true;  // catches also pyridazines, 1,2,3-triazines, etc.
              fg[fg_oxohetarene] := true;
              if opt_pos then add2fgloc(fg_oxohetarene,a_view);   // v0.5
            end;
          if (is_amino(a_view,a_n)) or
             ((atom^[a_n].atype = 'NAM') and (atom^[a_n].neighbor_count = 1) ) then 
            begin
              fg[fg_carboxylic_acid_amide]      := true;
              if opt_pos then add2fgloc(fg_carboxylic_acid_amide,a_view);   // v0.5
              fg[fg_carboxylic_acid_prim_amide] := true;
              if opt_pos then add2fgloc(fg_carboxylic_acid_prim_amide,a_view);   // v0.5
            end;
          //if (is_alkylamino(a_view,a_n)) or (is_arylamino(a_view,a_n)) then 
          if (is_C_monosubst_amino(a_view,a_n)) and (not is_subst_acylamino(a_view,a_n)) then   // v0.3j
            begin
              if (bond^[get_bond(a_view,a_n)].arom = false) then 
                begin
                  fg[fg_carboxylic_acid_amide]      := true;
                  if opt_pos then add2fgloc(fg_carboxylic_acid_amide,a_view);   // v0.5
                end;
              if (bond^[get_bond(a_view,a_n)].arom = false) then 
                begin
                  fg[fg_carboxylic_acid_sec_amide]  := true;
                  if opt_pos then add2fgloc(fg_carboxylic_acid_sec_amide,a_view);   // v0.5
                end;
              if (bond^[get_bond(a_view,a_n)].ring_count > 0) then
                begin
                  if (bond^[get_bond(a_view,a_n)].arom = true) then
                    begin
                     fg[fg_oxohetarene]    := true;
                     if opt_pos then add2fgloc(fg_oxohetarene,a_view);   // v0.5
                    end else 
                      begin
                        fg[fg_lactam]               := true;
                        if opt_pos then add2fgloc(fg_lactam,a_view);   // v0.5
                      end;
                end;
            end;          
          //if (is_dialkylamino(a_view,a_n)) or (is_alkylarylamino(a_view,a_n)) or
          //   (is_diarylamino(a_view,a_n)) then 
          if (is_C_disubst_amino(a_view,a_n)) and (not is_subst_acylamino(a_view,a_n)) then   // v0.3j
            begin
              if (bond^[get_bond(a_view,a_n)].arom = false) then 
                begin
                  fg[fg_carboxylic_acid_amide]      := true;
                  if opt_pos then add2fgloc(fg_carboxylic_acid_amide,a_view);   // v0.5
                end;
              if (bond^[get_bond(a_view,a_n)].arom = false) then 
                begin
                  fg[fg_carboxylic_acid_tert_amide] := true;
                  if opt_pos then add2fgloc(fg_carboxylic_acid_tert_amide,a_view);   // v0.5
                end;
              if (bond^[get_bond(a_view,a_n)].ring_count > 0) then
                begin
                  if (bond^[get_bond(a_view,a_n)].arom = true) then
                    //fg[fg_lactam_heteroarom]    := true else 
                    begin
                      fg[fg_oxohetarene]    := true;
                      if opt_pos then add2fgloc(fg_oxohetarene,a_view);   // v0.5
                    end else 
                      begin
                        fg[fg_lactam]               := true;
                        if opt_pos then add2fgloc(fg_lactam,a_view);   // v0.5
                      end;
                end;
            end;
          if (is_hydroxylamino(a_view,a_n)) then 
            begin
              fg[fg_hydroxamic_acid]            := true;
              if opt_pos then add2fgloc(fg_hydroxamic_acid,a_view);   // v0.5
            end;
          if (is_hydrazino(a_view,a_n)) then 
            begin
              fg[fg_carboxylic_acid_hydrazide]  := true;
              if opt_pos then add2fgloc(fg_carboxylic_acid_hydrazide,a_view);   // v0.5
            end;
          if (is_azido(a_view,a_n)) then 
            begin
              fg[fg_carboxylic_acid_azide]  := true;
              if opt_pos then add2fgloc(fg_carboxylic_acid_azide,a_view);   // v0.5
            end;
        end;
      if (s_count = 1) then 
        begin  // anhydride is checked somewhere else
          if (bond^[get_bond(a_view,a_s)].arom = false) then 
            begin
              fg[fg_thiocarboxylic_acid_deriv]  := true;
              if opt_pos then add2fgloc(fg_thiocarboxylic_acid_deriv,a_view);   // v0.5
            end;
          if (is_sulfanyl(a_view,a_s)) then 
            begin
              fg[fg_thiocarboxylic_acid] := true;
              if opt_pos then add2fgloc(fg_thiocarboxylic_acid,a_view);   // v0.5
            end;
          if (is_alkylsulfanyl(a_view,a_s)) or (is_arylsulfanyl(a_view,a_s)) then 
            begin
              if (bond^[get_bond(a_view,a_s)].arom = false) then 
                begin
                  fg[fg_thiocarboxylic_acid_ester] := true;
                  if opt_pos then add2fgloc(fg_thiocarboxylic_acid_ester,a_view);   // v0.5
                end;
              if (bond^[get_bond(a_view,a_s)].ring_count > 0) then
                begin
                  if (bond^[get_bond(a_view,a_s)].arom = true) then
                    //fg[fg_thiolactone_heteroarom] := true else fg[fg_thiolactone] := true;
                    begin
                      fg[fg_oxohetarene] := true;
                      if opt_pos then add2fgloc(fg_oxohetarene,a_view);   // v0.5
                    end else 
                    begin
                      fg[fg_thiolactone] := true;
                      if opt_pos then add2fgloc(fg_thiolactone,a_view);   // v0.5
                    end;
                end;
            end;            
        end;
    end;  // end oxo-C
  if (is_thioxo_c(a_view)) then
    begin
      // fg[fg_thiocarboxylic_acid_deriv]  := true;
      if (o_count = 1) then 
        begin  // anhydride is checked somewhere else
          if (bond^[get_bond(a_view,a_o)].arom = false) then 
            begin
              fg[fg_thiocarboxylic_acid_deriv]  := true;
              if opt_pos then add2fgloc(fg_thiocarboxylic_acid_deriv,a_view);   // v0.5
            end;
          if (is_hydroxy(a_view,a_o)) then 
            begin
              fg[fg_thiocarboxylic_acid] := true;   // fixed in v0.3c
              if opt_pos then add2fgloc(fg_thiocarboxylic_acid,a_view);   // v0.5
            end;
          if (is_alkoxy(a_view,a_o)) or (is_aryloxy(a_view,a_o)) then 
            begin
              if (bond^[get_bond(a_view,a_s)].arom = false) then 
                begin
                  fg[fg_thiocarboxylic_acid_ester] := true;
                  if opt_pos then add2fgloc(fg_thiocarboxylic_acid_ester,a_view);   // v0.5
                end;
              if (bond^[get_bond(a_view,a_o)].ring_count > 0) then
                begin
                  if (bond^[get_bond(a_view,a_o)].arom = true) then
                    //fg[fg_thiolactone_heteroarom] := true else fg[fg_thiolactone] := true;
                    begin
                      fg[fg_thioxohetarene] := true;
                      if opt_pos then add2fgloc(fg_thioxohetarene,a_view);   // v0.5
                    end else 
                    begin
                      fg[fg_thiolactone] := true;
                      if opt_pos then add2fgloc(fg_thiolactone,a_view);   // v0.5
                    end;
                end;
            end;            
        end;
      if (n_count = 1) then 
        begin
          if (bond^[get_bond(a_view,a_n)].arom = false) then
            begin
              fg[fg_thiocarboxylic_acid_deriv]  := true;
              if opt_pos then add2fgloc(fg_thiocarboxylic_acid_deriv,a_view);   // v0.5
            end else
            //fg[fg_thiolactam_heteroarom] := true;  // catches also pyridazines, 1,2,3-triazines, etc.
              begin  
                fg[fg_thioxohetarene] := true;  // catches also pyridazines, 1,2,3-triazines, etc.
                if opt_pos then add2fgloc(fg_thioxohetarene,a_view);   // v0.5
              end;
          if ((is_amino(a_view,a_n)) or (atom^[a_n].neighbor_count = 1)) then   // v0.4a
            begin
              fg[fg_thiocarboxylic_acid_amide]      := true;
              if opt_pos then add2fgloc(fg_thiocarboxylic_acid_amide,a_view);   // v0.5
              // fg[fg_thiocarboxylic_acid_prim_amide] := true;
            end;
          //if (is_alkylamino(a_view,a_n)) or (is_arylamino(a_view,a_n)) then 
          if (is_C_monosubst_amino(a_view,a_n)) and (not is_subst_acylamino(a_view,a_n)) then   // v0.3j
            begin
              if (bond^[get_bond(a_view,a_n)].arom = false) then 
                begin
                  fg[fg_thiocarboxylic_acid_amide]      := true;
                  if opt_pos then add2fgloc(fg_thiocarboxylic_acid_amide,a_view);   // v0.5
                end;
              //fg[fg_thiocarboxylic_acid_sec_amide]  := true;
              if (bond^[get_bond(a_view,a_n)].ring_count > 0) then
                begin
                  if (bond^[get_bond(a_view,a_n)].arom = true) then
                    //fg[fg_thiolactam_heteroarom] := true else fg[fg_thiolactam] := true;
                    begin
                      fg[fg_thioxohetarene] := true;
                      if opt_pos then add2fgloc(fg_thioxohetarene,a_view);   // v0.5
                    end else 
                    begin
                      fg[fg_thiolactam] := true;
                      if opt_pos then add2fgloc(fg_thiolactam,a_view);   // v0.5
                    end;
                end;
            end;          
          //if (is_dialkylamino(a_view,a_n)) or (is_alkylarylamino(a_view,a_n)) or
          //   (is_diarylamino(a_view,a_n)) then 
          if (is_C_disubst_amino(a_view,a_n)) and (not is_subst_acylamino(a_view,a_n)) then   // v0.3j
            begin
              if (bond^[get_bond(a_view,a_n)].arom = false) then 
                begin
                  fg[fg_thiocarboxylic_acid_amide]      := true;
                  if opt_pos then add2fgloc(fg_thiocarboxylic_acid_amide,a_view);   // v0.5
                end;
              //fg[fg_thiocarboxylic_acid_tert_amide] := true;
              if (bond^[get_bond(a_view,a_n)].ring_count > 0) then
                begin
                  if (bond^[get_bond(a_view,a_n)].arom = true) then
                    //fg[fg_thiolactam_heteroarom] := true else fg[fg_thiolactam] := true;
                    begin
                      fg[fg_thioxohetarene] := true;
                      if opt_pos then add2fgloc(fg_thioxohetarene,a_view);   // v0.5
                    end else 
                    begin
                      fg[fg_thiolactam] := true;
                      if opt_pos then add2fgloc(fg_thiolactam,a_view);   // v0.5
                    end;
                end;
            end;          
        end;
      if (s_count = 1) then
        begin  // anhydride is checked somewhere else
          if (bond^[get_bond(a_view,a_s)].arom = false) then 
            begin
              fg[fg_thiocarboxylic_acid_deriv]  := true;
              if opt_pos then add2fgloc(fg_thiocarboxylic_acid_deriv,a_view);   // v0.5
            end;
          if (is_sulfanyl(a_view,a_s)) then 
            begin
              fg[fg_thiocarboxylic_acid] := true;
              if opt_pos then add2fgloc(fg_thiocarboxylic_acid,a_view);   // v0.5
            end;
          if (is_alkylsulfanyl(a_view,a_s)) or (is_arylsulfanyl(a_view,a_s)) then 
            begin
              if (bond^[get_bond(a_view,a_s)].arom = false) then 
                begin
                  fg[fg_thiocarboxylic_acid_ester] := true;
                  if opt_pos then add2fgloc(fg_thiocarboxylic_acid_ester,a_view);   // v0.5
                end;
              if (bond^[get_bond(a_view,a_s)].ring_count > 0) then
                begin
                  if (bond^[get_bond(a_view,a_s)].arom = true) then
                    //fg[fg_thiolactone_heteroarom] := true else fg[fg_thiolactone] := true;
                    begin
                      fg[fg_thioxohetarene] := true;
                      if opt_pos then add2fgloc(fg_thioxohetarene,a_view);   // v0.5
                    end else 
                    begin
                      fg[fg_thiolactone] := true;
                      if opt_pos then add2fgloc(fg_thiolactone,a_view);   // v0.5
                    end;
                end;
            end;            
        end;
    end;  // end thioxo-C
  if (is_true_imino_c(a_view)) then
    begin
      if (o_count = 1) then 
        begin  
          if (bond^[get_bond(a_view,a_o)].arom = false) then 
            begin
              fg[fg_carboxylic_acid_deriv]  := true;
              if opt_pos then add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
            end;
          if (is_alkoxy(a_view,a_o)) or (is_aryloxy(a_view,a_o)) then 
            begin
              if (bond^[get_bond(a_view,a_o)].arom = false) then
                begin
                  fg[fg_imido_ester] := true;
                  if opt_pos then add2fgloc(fg_imido_ester,a_view);   // v0.5
                end;
            end;            
        end;
      if (n_count = 1) and (bond^[get_bond(a_view,a_n)].arom = false) then 
        begin
          if (bond^[get_bond(a_view,a_n)].arom = false) then 
            begin
              fg[fg_carboxylic_acid_deriv]  := true;
              if opt_pos then add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
            end;
          if (is_amino(a_view,a_n)) or (is_subst_amino(a_view,a_n)) then
            begin
              if (bond^[get_bond(a_view,a_n)].arom = false) then 
                begin
                  fg[fg_carboxylic_acid_deriv]  := true;
                  if opt_pos then add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
                  fg[fg_carboxylic_acid_amidine]      := true;
                  if opt_pos then add2fgloc(fg_carboxylic_acid_amidine,a_view);   // v0.5
                end;
            end;          
          if (is_hydrazino(a_view,a_n)) then 
            begin
              if (bond^[get_bond(a_view,a_n)].arom = false) then 
                begin
                  fg[fg_carboxylic_acid_amidrazone]  := true;
                  if opt_pos then 
                    begin
                      add2fgloc(fg_carboxylic_acid_amidrazone,a_view);   // v0.5
                      fgloc_set_hydrazino(a_n);   // v0.5
                    end;
                end;
            end;
        end;
      if (n_count = 1) and (bond^[get_bond(a_view,a_n)].arom = true) then // catches also pyridazines, 1,2,3-triazines, etc.
        begin
          fg[fg_iminohetarene] := true;
          if opt_pos then add2fgloc(fg_iminohetarene,a_view);   // v0.5
        end;
      if (s_count = 1) then 
        begin  
          if (bond^[get_bond(a_view,a_s)].arom = false) then 
            begin
              fg[fg_carboxylic_acid_deriv]  := true;
              if opt_pos then add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
            end;
          if (is_alkylsulfanyl(a_view,a_s)) or (is_arylsulfanyl(a_view,a_s)) then 
            begin
              if (bond^[get_bond(a_view,a_s)].arom = false) then
                begin
                  fg[fg_imido_thioester] := true;
                  if opt_pos then add2fgloc(fg_imido_thioester,a_view);   // v0.5
                end;
            end;            
        end;
    end;   // is true_imino_c
  if is_hydroximino_c(a_view) then
    begin
      if (bond^[get_bond(a_view,a_n)].arom = false) then 
        begin
          fg[fg_carboxylic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
        end;
      if (o_count = 1) then 
        begin  
          if (is_hydroxy(a_view,a_o)) then 
            begin
              fg[fg_hydroxamic_acid] := true;
              if opt_pos then add2fgloc(fg_hydroxamic_acid,a_view);   // v0.5
            end;            
        end;
    end;   // is hydroximino_c
  if is_hydrazono_c(a_view) then
    begin
      if (bond^[get_bond(a_view,a_n)].arom = false) then 
        begin
          fg[fg_carboxylic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_carboxylic_acid_deriv,a_view);   // v0.5
        end;
      if (n_count = 1) then 
        begin  
          if (is_amino(a_view,a_n)) or
             (is_subst_amino(a_view,a_n)) then 
            begin
              fg[fg_carboxylic_acid_amidrazone] := true;
              if opt_pos then add2fgloc(fg_carboxylic_acid_amidrazone,a_view);   // v0.5
            end;            
        end;
    end;   // is hydrazono_c
end;


procedure chk_co2_sp2(a_view,a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  o_count, or_count, n_count, nn_count, nnx_count, s_count, sr_count : integer;
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_view);
  o_count := 0; or_count := 0; n_count := 0;
  nn_count := 0; nnx_count := 0; s_count := 0; sr_count := 0;
  for i := 1 to atom^[a_view].neighbor_count do
    begin
      if bond^[(get_bond(a_view,nb[i]))].btype = 'S' then
        begin
          nb_el := atom^[(nb[i])].element;
          if (nb_el <> 'C ') then
            begin
              if (nb_el = 'O ') then
                begin
                  inc(o_count);
                  if (is_alkoxy(a_view,nb[i])) or (is_alkenyloxy(a_view,nb[i])) or  // v0.3j
                     (is_aryloxy(a_view,nb[i])) then inc (or_count);
                end;
              if (nb_el = 'N ') then 
                begin 
                  inc(n_count); 
                  if (is_hydrazino(a_view,nb[i])) then inc(nn_count);
                  if (is_subst_hydrazino(a_view,nb[i])) then inc(nnx_count);  // more general...
                end;
              if (nb_el = 'S ') then 
                begin 
                  inc(s_count); 
                  if (is_alkylsulfanyl(a_view,nb[i])) or 
                     (is_arylsulfanyl(a_view,nb[i])) then inc (sr_count);
                end;
            end;
        end;
    end;     
  if (is_oxo_c(a_view)) then
    begin
      if (o_count = 2) then 
        begin  
          fg[fg_carbonic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_carbonic_acid_deriv,a_view);   // v0.5
          if (or_count = 1) then 
            begin
              fg[fg_carbonic_acid_monoester]  := true;
              if opt_pos then add2fgloc(fg_carbonic_acid_monoester,a_view);   // v0.5
            end;
          if (or_count = 2) then 
            begin
              fg[fg_carbonic_acid_diester]  := true;
              if opt_pos then add2fgloc(fg_carbonic_acid_diester,a_view);   // v0.5
            end;
        end;
      if (o_count = 1) and (s_count = 1) then 
        begin  
          fg[fg_thiocarbonic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_thiocarbonic_acid_deriv,a_view);   // v0.5
          if (or_count + sr_count = 1) then 
            begin
              fg[fg_thiocarbonic_acid_monoester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbonic_acid_monoester,a_view);   // v0.5
            end;
          if (or_count + sr_count = 2) then 
            begin
              fg[fg_thiocarbonic_acid_diester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbonic_acid_diester,a_view);   // v0.5
            end;
        end;
      if (s_count = 2) then 
        begin  
          fg[fg_thiocarbonic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_thiocarbonic_acid_deriv,a_view);   // v0.5
          if (sr_count = 1) then 
            begin
              fg[fg_thiocarbonic_acid_monoester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbonic_acid_monoester,a_view);   // v0.5
            end;
          if (sr_count = 2) then 
            begin
              fg[fg_thiocarbonic_acid_diester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbonic_acid_diester,a_view);   // v0.5
            end;
        end;
      if (o_count = 1) and (n_count = 1) then 
        begin  
          fg[fg_carbamic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_carbamic_acid_deriv,a_view);   // v0.5
          if (or_count = 0) then 
            begin
              fg[fg_carbamic_acid]  := true;
              if opt_pos then add2fgloc(fg_carbamic_acid,a_view);   // v0.5
            end;
          if (or_count = 1) then 
            begin
              fg[fg_carbamic_acid_ester]  := true;
              if opt_pos then add2fgloc(fg_carbamic_acid_ester,a_view);   // v0.5
            end;
        end;
      if (s_count = 1) and (n_count = 1) then 
        begin  
          fg[fg_thiocarbamic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_thiocarbamic_acid_deriv,a_view);   // v0.5
          if (sr_count = 0) then 
            begin
              fg[fg_thiocarbamic_acid]  := true;
              if opt_pos then add2fgloc(fg_thiocarbamic_acid,a_view);   // v0.5
            end;
          if (sr_count = 1) then 
            begin
              fg[fg_thiocarbamic_acid_ester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbamic_acid_ester,a_view);   // v0.5
            end;
        end;
      if (n_count = 2) then
        begin  
          if (nn_count = 1) then 
            begin
              fg[fg_semicarbazide]  := true;
              if opt_pos then add2fgloc(fg_semicarbazide,a_view);   // v0.5
            end else
            begin
              if (nnx_count = 0) then 
                begin
                  fg[fg_urea] := true;  // excludes semicarbazones
                  if opt_pos then add2fgloc(fg_urea,a_view);   // v0.5
                end;
            end;                     
        end;  
    end;  // end oxo-C
  if (is_thioxo_c(a_view)) then
    begin
      if (o_count = 2) then 
        begin  
          fg[fg_thiocarbonic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_thiocarbonic_acid_deriv,a_view);   // v0.5
          if (or_count = 1) then 
            begin
              fg[fg_thiocarbonic_acid_monoester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbonic_acid_monoester,a_view);   // v0.5
            end;
          if (or_count = 2) then 
            begin
              fg[fg_thiocarbonic_acid_diester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbonic_acid_diester,a_view);   // v0.5
            end;
        end;
      if (o_count = 1) and (s_count = 1) then 
        begin  
          fg[fg_thiocarbonic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_thiocarbonic_acid_deriv,a_view);   // v0.5
          if (or_count + sr_count = 1) then 
            begin
              fg[fg_thiocarbonic_acid_monoester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbonic_acid_monoester,a_view);   // v0.5
            end;
          if (or_count + sr_count = 2) then 
            begin
              fg[fg_thiocarbonic_acid_diester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbonic_acid_diester,a_view);   // v0.5
            end;
        end;
      if (s_count = 2) then 
        begin  
          fg[fg_thiocarbonic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_thiocarbonic_acid_deriv,a_view);   // v0.5
          if (sr_count = 1) then 
            begin
              fg[fg_thiocarbonic_acid_monoester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbonic_acid_monoester,a_view);   // v0.5
            end;
          if (sr_count = 2) then 
            begin
              fg[fg_thiocarbonic_acid_diester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbonic_acid_diester,a_view);   // v0.5
            end;
        end;
      if (o_count = 1) and (n_count = 1) then 
        begin  
          fg[fg_thiocarbamic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_thiocarbamic_acid_deriv,a_view);   // v0.5
          if (or_count = 0) then 
            begin
              fg[fg_thiocarbamic_acid]  := true;
              if opt_pos then add2fgloc(fg_thiocarbamic_acid,a_view);   // v0.5
            end;
          if (or_count = 1) then 
            begin
              fg[fg_thiocarbamic_acid_ester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbamic_acid_ester,a_view);   // v0.5
            end;
        end;
      if (s_count = 1) and (n_count = 1) then 
        begin  
          fg[fg_thiocarbamic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_thiocarbamic_acid_deriv,a_view);   // v0.5
          if (sr_count = 0) then 
            begin
              fg[fg_thiocarbamic_acid]  := true;
              if opt_pos then add2fgloc(fg_thiocarbamic_acid,a_view);   // v0.5
            end;
          if (sr_count = 1) then 
            begin
              fg[fg_thiocarbamic_acid_ester]  := true;
              if opt_pos then add2fgloc(fg_thiocarbamic_acid_ester,a_view);   // v0.5
            end;
        end;
      if (n_count = 2) then 
        begin  
          if (nn_count = 1) then 
            begin
              fg[fg_thiosemicarbazide]  := true;
              if opt_pos then add2fgloc(fg_thiosemicarbazide,a_view);   // v0.5
            end else
            begin
              if (nnx_count = 0) then 
                begin
                  fg[fg_thiourea] := true;  // excludes thiosemicarbazones
                  if opt_pos then add2fgloc(fg_thiourea,a_view);   // v0.5
                end;
            end;                     
        end;  
    end;  // end thioxo-C
  if (is_true_imino_c(a_view)) and (bond^[get_bond(a_view,a_ref)].arom = false) then
    begin
      if (o_count = 1) and (n_count = 1) then 
        begin  
          fg[fg_isourea]  := true;
          if opt_pos then add2fgloc(fg_isourea,a_view);   // v0.5
        end;
      if (s_count = 1) and (n_count = 1) then 
        begin  
          fg[fg_isothiourea]  := true;
          if opt_pos then add2fgloc(fg_isothiourea,a_view);   // v0.5
        end;
      if (n_count = 2) then 
        begin  
          fg[fg_guanidine] := true;
          if opt_pos then add2fgloc(fg_guanidine,a_view);   // v0.5
        end;  
    end;  // end Imino-C
end;


procedure chk_co2_sp(a_view,a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  o_count, n_count, s_count : integer;
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_view);
  o_count := 0; n_count := 0; s_count := 0;
  for i := 1 to atom^[a_view].neighbor_count do
    begin
      if bond^[(get_bond(a_view,nb[i]))].btype = 'D' then
        begin
          nb_el := atom^[(nb[i])].element;
          if (nb_el <> 'C ') then 
            begin
              if (nb_el = 'O ') then 
                begin 
                  inc(o_count); 
                end;
              if (nb_el = 'N ') then 
                begin 
                  inc(n_count); 
                end;
              if (nb_el = 'S ') then 
                begin 
                  inc(s_count); 
                end;
            end;
        end;
    end;     
  if (o_count + s_count = 2) then 
    begin
      fg[fg_co2_deriv] := true;  // new in v0.3b
      if opt_pos then add2fgloc(fg_co2_deriv,a_view);   // v0.5
    end;
  if (o_count = 1) and (n_count = 1) then 
    begin
      fg[fg_isocyanate] := true;
      if opt_pos then add2fgloc(fg_isocyanate,a_view);   // v0.5
    end;
  if (s_count = 1) and (n_count = 1) then 
    begin
      fg[fg_isothiocyanate] := true;
      if opt_pos then add2fgloc(fg_isothiocyanate,a_view);   // v0.5
    end;
  if (n_count = 2) then 
    begin
      fg[fg_carbodiimide] := true;
      if opt_pos then add2fgloc(fg_carbodiimide,a_view);   // v0.5
    end;
end;


procedure chk_triple(a1,a2:integer);
var
  a1_el, a2_el : str2;
  b_id : integer;  // v0.5
begin
  a1_el := atom^[a1].element;
  a2_el := atom^[a2].element;
  if (a1_el = 'C ') and (a2_el = 'C ') and (bond^[get_bond(a1,a2)].arom = false) then
    begin
      fg[fg_alkyne] := true;
      if opt_pos then
        begin
          b_id := get_bond(a1,a2);
          add2fgloc(fg_alkyne,b_id);
        end;
    end;
  if (is_nitrile(a1,a2))     then 
    begin
      fg[fg_nitrile]     := true;
      if opt_pos then add2fgloc(fg_nitrile,a1);
    end;
  if (is_isonitrile(a1,a2))  then 
    begin
      fg[fg_isonitrile]  := true;
      if opt_pos then add2fgloc(fg_isonitrile,a2);  // a2 is the N atom!
    end;
  if (is_cyanate(a1,a2))     then 
    begin
      fg[fg_cyanate]     := true;
      if opt_pos then add2fgloc(fg_cyanate,a1);
    end;
  if (is_thiocyanate(a1,a2)) then 
    begin
      fg[fg_thiocyanate] := true;
      if opt_pos then add2fgloc(fg_thiocyanate,a1);
    end;
end;


procedure chk_ccx(a_view,a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  oh_count, or_count, n_count : integer;
  b_id : integer;  // v0.5
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_ref);
  oh_count := 0; or_count := 0; n_count := 0;
  for i := 1 to atom^[a_ref].neighbor_count do
    begin
      if bond^[(get_bond(a_ref,nb[i]))].btype = 'S' then
        begin
          if (is_hydroxy(a_ref,nb[i])) then inc(oh_count);
          if (is_alkoxy(a_ref,nb[i])) or 
             (is_aryloxy(a_ref,nb[i])) or
             (is_siloxy(a_ref,nb[i])) then inc(or_count);
          if (atom^[(nb[i])].atype = 'N3 ') or 
             (atom^[(nb[i])].atype = 'NAM') then inc(n_count);
        end;
    end;
  if (oh_count = 1) then 
    begin
      fg[fg_enol]      := true;
      if opt_pos then add2fgloc(fg_enol,a_ref);  // v0.5
    end;
  if (or_count = 1) then 
    begin
      fg[fg_enolether] := true;  
      if opt_pos then add2fgloc(fg_enolether,a_ref);  // v0.5
    end;
  if (n_count = 1)  then 
    begin
      fg[fg_enamine]   := true;
      if opt_pos then add2fgloc(fg_enamine,a_ref);  // v0.5
    end;
  // new in v0.2f   (regard anything else as an alkene)
  if ((oh_count + or_count + n_count) = 0) then 
    begin
      fg[fg_alkene] := true;
        if opt_pos then   // v0.5
          begin
            b_id := get_bond(a_view,a_ref);
            add2fgloc(fg_alkene,b_id);
          end;
    end;
end;


procedure chk_xccx(a_view,a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  oh_count, or_count, n_count : integer;
  b_id : integer;  // v0.5
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_view);
  oh_count := 0; or_count := 0; n_count := 0;
  for i := 1 to atom^[a_view].neighbor_count do
    begin
      if bond^[(get_bond(a_view,nb[i]))].btype = 'S' then
        begin
          if (is_hydroxy(a_view,nb[i])) then inc(oh_count);
          if (is_alkoxy(a_view,nb[i])) or
             (is_aryloxy(a_view,nb[i])) or
             (is_siloxy(a_view,nb[i])) then inc(or_count);
          if (atom^[(nb[i])].atype = 'N3 ') or
             (atom^[(nb[i])].atype = 'NAM') then inc(n_count);
        end;
    end;
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_ref);
  for i := 1 to atom^[a_ref].neighbor_count do
    begin
      if bond^[(get_bond(a_ref,nb[i]))].btype = 'S' then
        begin
          if (is_hydroxy(a_ref,nb[i])) then inc(oh_count);
          if (is_alkoxy(a_ref,nb[i])) or
             (is_aryloxy(a_ref,nb[i])) or
             (is_siloxy(a_ref,nb[i])) then inc(or_count);
          if (atom^[(nb[i])].atype = 'N3 ') or
             (atom^[(nb[i])].atype = 'NAM') then inc(n_count);
        end;
    end;
  if (oh_count = 2) then 
    begin
      fg[fg_enediol]      := true;
      if opt_pos then   // v0.5
        begin
          b_id := get_bond(a_view,a_ref);
          add2fgloc(fg_enediol,b_id);
        end;
    end;
  // new in v0.2f   (regard anything else as an alkene)
  if ((oh_count + or_count + n_count) = 0) then 
    begin
      fg[fg_alkene] := true;
      if opt_pos then   // v0.5
        begin
          b_id := get_bond(a_view,a_ref);
          add2fgloc(fg_alkene,b_id);
        end;
    end;
end;


procedure chk_n_o_dbl(a1,a2:integer);
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  or_count, n_count, c_count : integer;
  b : integer;          // v0.3j
  het_count : integer;  // v0.3k
  bt : char;            // v0.3k
  bo_sum : single;      // v0.3k
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a1);
  or_count := 0; n_count := 0; c_count := 0;
  bo_sum := 0;    // v0.3k
  het_count := 0; // v0.3k
  for i := 1 to atom^[a1].neighbor_count do
    begin
      if (nb[i] <> a2) then
        begin
          b := get_bond(a1,nb[i]);         // v0.3j
          nb_el := atom^[(nb[i])].element;
          bt := bond^[b].btype;                           // v0.3k
          if ((nb_el <> 'C ') and (nb_el <> 'H ') and (nb_el <> 'D ') and  // added 'D ' in v0.3n
             (nb_el <> 'DU') and (nb_el <> 'LP')) and     // v0.3k: ignore hetero atoms
             (bond^[b].arom = false) then inc(het_count); // in aromatic rings like isoxazole 
          if bt = 'S' then bo_sum := bo_sum + 1;           
          if bt = 'D' then bo_sum := bo_sum + 2;
          if bt = 'A' then bo_sum := bo_sum + 1.5;
          if (nb_el = 'O ') then inc(or_count);
          if (nb_el = 'N ') then inc(n_count);
          if (nb_el = 'C ') and (bond^[b].btype = 'S') then inc(c_count);  // v0.3k
          // if (is_alkyl(a1,nb[i])) or (is_aryl(a1,nb[i])) then inc(c_count);
        end;
    end;
  if ((or_count + n_count + c_count) = 1) and (atom^[a1].neighbor_count = 2) then  // excludes nitro etc.
    begin
      if (or_count = 1) then 
        begin
          fg[fg_nitrite]       := true;
          if opt_pos then add2fgloc(fg_nitrite,a1);   // v0.5
        end;
      if (c_count = 1)  then 
        begin
          fg[fg_nitroso_compound] := true;
          if opt_pos then add2fgloc(fg_nitroso_compound,a1);   // v0.5
        end;
      if (n_count = 1) then 
        begin
          fg[fg_nitroso_compound]  := true; // instead of nitrosamine  v0.3j
          if opt_pos then add2fgloc(fg_nitroso_compound,a1);   // v0.5
        end;
      //if (n_count = 1) then fg[fg_nitrosamine]   := true;  // still missing
    end;
  //if ((c_count > 1) and (or_count = 0) and (n_count = 0)) then
  //  begin
  //    fg[fg_n_oxide] := true;
  //  end;
  // new approach in v0.3k
  if ((het_count = 0) and (bo_sum > 2)) then   // =O does not count!
    begin
      fg[fg_n_oxide] := true;
      if opt_pos then add2fgloc(fg_n_oxide,a1);   // v0.5
    end;
end;


procedure chk_sulfoxide(a1,a2:integer);
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  o_count, c_count : integer;
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a1);
  o_count := 0; c_count := 0;
  for i := 1 to atom^[a1].neighbor_count do
    begin
      nb_el := atom^[(nb[i])].element;
      if (nb_el = 'O ') then inc(o_count);
      if (is_alkyl(a1,nb[i]) or is_aryl(a1,nb[i]) or 
         is_alkenyl(a1,nb[i]) or is_alkynyl(a1,nb[i])) then inc(c_count);   // v0.5
    end;
  if (o_count = 1) and (c_count = 2) then 
    begin
      fg[fg_sulfoxide] := true;
      if opt_pos then add2fgloc(fg_sulfoxide,a1);   // v0.5
    end;
end;


procedure chk_double(a1,a2:integer);
var
  a1_el, a2_el : str2;
  b_id : integer;  // v0.5
begin
  a1_el := atom^[a1].element;
  a2_el := atom^[a2].element;
  if (a1_el = 'C ') and (a2_el <> 'C ') and (bond^[get_bond(a1,a2)].arom = false) then
    begin
      if (hetbond_count(a1) = 2) then
        begin
          chk_carbonyl_deriv(a1,a2);
        end;
      if (hetbond_count(a1) = 3) then
        begin
          chk_carboxyl_deriv(a1,a2);
        end;
      if (hetbond_count(a1) = 4) then
        begin
          if (atom^[a1].atype = 'C2 ') then chk_co2_sp2(a1,a2);
          if (atom^[a1].atype = 'C1 ') then chk_co2_sp(a1,a2);
        end;
    end;  // end C=X
  if (atom^[a1].atype = 'C2 ') and (atom^[a2].atype = 'C2 ') and 
     (bond^[get_bond(a1,a2)].arom = false) then
    begin
      if (hetbond_count(a1) = 0) and (hetbond_count(a2) = 2) then 
        begin
          fg[fg_ketene_acetal_deriv] := true;
          if opt_pos then add2fgloc(fg_ketene_acetal_deriv,a2);  // v0.5
        end;
      if (hetbond_count(a1) = 0) and (hetbond_count(a2) = 1) then chk_ccx(a1,a2);
      if (hetbond_count(a1) = 1) and (hetbond_count(a2) = 1) then chk_xccx(a1,a2);
      if (hetbond_count(a1) = 0) and (hetbond_count(a2) = 0) and
         (atom^[a1].arom = false) and (atom^[a2].arom = false) then 
        begin
          fg[fg_alkene] := true;
          if opt_pos then   // v0.5
            begin
              b_id := get_bond(a1,a2);
              add2fgloc(fg_alkene,b_id);
            end;
        end;
    end;
  if (a1_el = 'N ') and (a2_el = 'N ') and 
     (hetbond_count(a1) = 2) and (hetbond_count(a2) = 2) and
     (bond^[get_bond(a1,a2)].arom = false) and
     (atom^[a1].neighbor_count = 2) and (atom^[a2].neighbor_count = 2) 
       then 
         begin
           fg[fg_azo_compound] := true;
           if opt_pos then  // v0.5
             begin
               b_id := get_bond(a1,a2);
               add2fgloc(fg_azo_compound,b_id);
             end;
         end;
  if (a1_el = 'N ') and (a2_el = 'O ') then chk_n_o_dbl(a1,a2);
  if (a1_el = 'S ') and (a2_el = 'O ') then chk_sulfoxide(a1,a2);
end;


procedure chk_c_hal(a1,a2:integer);
var
  a2_el : str2;
begin
  a2_el := atom^[a2].element;
  fg[fg_halogen_deriv] := true;
  if opt_pos then add2fgloc(fg_halogen_deriv,a2);  // v0.5
  if atom^[a1].arom then begin
    fg[fg_aryl_halide] := true;
    if (a2_el = 'F ') then 
      begin
        fg[fg_aryl_fluoride] := true;
        if opt_pos then add2fgloc(fg_aryl_fluoride,a2);  // v0.5
      end;
    if (a2_el = 'CL') then 
      begin
        fg[fg_aryl_chloride] := true;
        if opt_pos then add2fgloc(fg_aryl_chloride,a2);  // v0.5
      end;
    if (a2_el = 'BR') then 
      begin
        fg[fg_aryl_bromide]  := true;
        if opt_pos then add2fgloc(fg_aryl_bromide,a2);  // v0.5
      end;
    if (a2_el = 'I ') then 
      begin
        fg[fg_aryl_iodide]   := true;
        if opt_pos then add2fgloc(fg_aryl_iodide,a2);  // v0.5
      end;
  end else
  begin
    if (atom^[a1].atype = 'C3 ') and (hetbond_count(a1) <= 2) then
      begin  // alkyl halides
        fg[fg_alkyl_halide] := true;
        if (a2_el = 'F ') then 
          begin
            fg[fg_alkyl_fluoride] := true;
            if opt_pos then add2fgloc(fg_alkyl_fluoride,a2);  // v0.5
          end;
        if (a2_el = 'CL') then 
          begin
            fg[fg_alkyl_chloride] := true;
            if opt_pos then add2fgloc(fg_alkyl_chloride,a2);  // v0.5
          end;
        if (a2_el = 'BR') then 
          begin
            fg[fg_alkyl_bromide]  := true;
            if opt_pos then add2fgloc(fg_alkyl_bromide,a2);  // v0.5
          end;
        if (a2_el = 'I ') then 
          begin
            fg[fg_alkyl_iodide]   := true;                      
            if opt_pos then add2fgloc(fg_alkyl_iodide,a2);  // v0.5
          end;
      end;
    if (atom^[a1].atype = 'C2 ') and (hetbond_count(a1) = 3) then
      begin  // acyl halides and related compounds
        if (is_oxo_c(a1)) then
          begin
            fg[fg_acyl_halide] := true;
            if opt_pos then add2fgloc(fg_acyl_halide,a1);  // v0.5
            if (a2_el = 'F ') then 
              begin
                fg[fg_acyl_fluoride] := true;
                if opt_pos then add2fgloc(fg_acyl_fluoride,a1);  // v0.5
              end;
            if (a2_el = 'CL') then 
              begin
                fg[fg_acyl_chloride] := true;
                if opt_pos then add2fgloc(fg_acyl_chloride,a1);  // v0.5
              end;
            if (a2_el = 'BR') then 
              begin
                fg[fg_acyl_bromide]  := true;
                if opt_pos then add2fgloc(fg_acyl_bromide,a1);  // v0.5
              end;
            if (a2_el = 'I ') then 
              begin
                fg[fg_acyl_iodide]   := true;                      
                if opt_pos then add2fgloc(fg_acyl_iodide,a1);  // v0.5
              end;
          end;
        if (is_thioxo_c(a1)) then
          begin
            fg[fg_thiocarboxylic_acid_deriv] := true;
            if opt_pos then add2fgloc(fg_thiocarboxylic_acid_deriv,a1);  // v0.5
          end;
        if (is_imino_c(a1)) then
          begin
            fg[fg_imidoyl_halide] := true;
            if opt_pos then add2fgloc(fg_imidoyl_halide,a1);  // v0.5
          end;
      end;
    if (atom^[a1].atype = 'C2 ') and (hetbond_count(a1) = 4) then
      begin  // chloroformates etc.
        fg[fg_co2_deriv] := true;
        if opt_pos then add2fgloc(fg_co2_deriv,a1);  // v0.5
        if (is_oxo_c(a1)) then
          begin
            fg[fg_carbonic_acid_deriv] := true;
            if opt_pos then add2fgloc(fg_carbonic_acid_deriv,a1);  // v0.5
            if (is_alkoxycarbonyl(a2,a1)) or (is_aryloxycarbonyl(a2,a1)) then 
              begin
                fg[fg_carbonic_acid_ester_halide] := true;
                if opt_pos then add2fgloc(fg_carbonic_acid_ester_halide,a1);  // v0.5
              end;
            if (is_carbamoyl(a2,a1)) then 
              begin 
                fg[fg_carbamic_acid_deriv]  := true;
                if opt_pos then add2fgloc(fg_carbamic_acid_deriv,a1);  // v0.5
                fg[fg_carbamic_acid_halide] := true;
                if opt_pos then add2fgloc(fg_carbamic_acid_halide,a1);  // v0.5
              end;
          end;
        if (is_thioxo_c(a1)) then
          begin
            fg[fg_thiocarbonic_acid_deriv] := true;
            if opt_pos then add2fgloc(fg_thiocarbonic_acid_deriv,a1);  // v0.5
            if (is_alkoxythiocarbonyl(a2,a1)) or (is_aryloxythiocarbonyl(a2,a1)) then
              begin 
                fg[fg_thiocarbonic_acid_ester_halide] := true;
                if opt_pos then add2fgloc(fg_thiocarbonic_acid_ester_halide,a1);  // v0.5
              end;
            if (is_thiocarbamoyl(a2,a1)) then 
              begin 
                fg[fg_thiocarbamic_acid_deriv]  := true;
                if opt_pos then add2fgloc(fg_thiocarbamic_acid_deriv,a1);  // v0.5
                fg[fg_thiocarbamic_acid_halide] := true;
                if opt_pos then add2fgloc(fg_thiocarbamic_acid_halide,a1);  // v0.5
              end;
          end;
      end;
    // still missing: polyhalogen compounds (-CX2H, -CX3)
  end;    // end of non-aromatic halogen compounds
end;


procedure chk_c_o(a1,a2:integer);
// a1 = C, a2 = O
begin
  // ignore heteroaromatic rings (like furan, thiophene, etc.)
  if (bond^[get_bond(a1,a2)].arom = true) then exit;
  if (is_true_alkyl(a2,a1)) and (is_hydroxy(a1,a2)) then
    begin
      fg[fg_hydroxy] := true;
      fg[fg_alcohol] := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_hydroxy,a2);
          add2fgloc(fg_alcohol,a2);
        end;
      if atom^[a1].neighbor_count <= 2 then 
        begin
          fg[fg_prim_alcohol] := true;
          if opt_pos then add2fgloc(fg_prim_alcohol,a2);   // v0.5
        end;
      if atom^[a1].neighbor_count = 3 then 
        begin
          fg[fg_sec_alcohol]  := true;
          if opt_pos then add2fgloc(fg_sec_alcohol,a2);   // v0.5
        end;
      if atom^[a1].neighbor_count = 4 then 
        begin
          fg[fg_tert_alcohol] := true;
          if opt_pos then add2fgloc(fg_tert_alcohol,a2);   // v0.5
        end;
    end;
  if (is_aryl(a2,a1)) and (is_hydroxy(a1,a2)) then
    begin
      fg[fg_hydroxy] := true;
      fg[fg_phenol]  := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_hydroxy,a2);
          add2fgloc(fg_phenol,a2);
        end;
    end;
  if (is_true_alkyl(a2,a1)) and (is_true_alkoxy(a1,a2)) then
    begin
      fg[fg_ether]        := true;
      fg[fg_dialkylether] := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_ether,a2);
          add2fgloc(fg_dialkylether,a2);
        end;
    end;
  if ((is_true_alkyl(a2,a1)) and (is_aryloxy(a1,a2))) or 
     ((is_aryl(a2,a1)) and (is_true_alkoxy(a1,a2))) then
    begin
      fg[fg_ether]        := true;
      fg[fg_alkylarylether] := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_ether,a2);
          add2fgloc(fg_alkylarylether,a2);
        end;
    end;
  if (is_aryl(a2,a1)) and (is_aryloxy(a1,a2)) then
    begin
      fg[fg_ether]        := true;
      fg[fg_diarylether] := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_ether,a2);
          add2fgloc(fg_diarylether,a2);
        end;
    end;
  if (is_true_alkyl(a2,a1) or is_aryl(a2,a1)) and (is_alkynyloxy(a1,a2)) then
    begin
      fg[fg_ether]        := true;
      ether_generic       := true;
      if opt_pos then add2fgloc(fg_ether,a2);   // v0.5
    end;
  if (is_alkynyl(a2,a1)) and (is_hydroxy(a1,a2)) then
    begin
      fg[fg_hydroxy]  := true;
      hydroxy_generic := true;
      if opt_pos then add2fgloc(fg_hydroxy,a2);   // v0.5
    end;
end;


procedure chk_c_s(a1,a2:integer);
// a1 = C, a2 = S
var
  i : integer;
  nb : neighbor_rec;
  nb_el : str2;
  o_count, oh_count, or_count, n_count, c_count, hal_count : integer;
begin
  // ignore heteroaromatic rings (like furan, thiophene, etc.)
  if (bond^[get_bond(a1,a2)].arom = true) then exit;
  if (is_alkyl(a2,a1)) and (is_sulfanyl(a1,a2)) then
    begin
      fg[fg_thiol] := true;
      fg[fg_alkylthiol] := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_thiol,a2);
          add2fgloc(fg_alkylthiol,a2);
        end;
    end;
  if (is_aryl(a2,a1)) and (is_sulfanyl(a1,a2)) then
    begin
      fg[fg_thiol]       := true;
      fg[fg_arylthiol]   := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_thiol,a2);
          add2fgloc(fg_arylthiol,a2);
        end;
    end;
  if (is_true_alkyl(a2,a1)) and (is_true_alkylsulfanyl(a1,a2)) then 
    begin
      fg[fg_thioether]  := true;
      if opt_pos then add2fgloc(fg_thioether,a2);   // v0.5
    end;
  if ((is_true_alkyl(a2,a1)) and (is_arylsulfanyl(a1,a2))) or 
     ((is_aryl(a2,a1)) and (is_true_alkylsulfanyl(a1,a2))) then 
    begin
      fg[fg_thioether] := true;
      if opt_pos then add2fgloc(fg_thioether,a2);   // v0.5
    end;
  if (is_aryl(a2,a1)) and (is_arylsulfanyl(a1,a2)) then 
    begin
      fg[fg_thioether]    := true;
      if opt_pos then add2fgloc(fg_thioether,a2);   // v0.5
    end;
  // check for sulfinic/sulfenic acid derivatives
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a2);
  o_count := 0; oh_count := 0; or_count := 0; n_count := 0; c_count := 0; hal_count := 0;
  for i := 1 to atom^[a2].neighbor_count do
    begin
      nb_el := atom^[(nb[i])].element;
      if (is_alkyl(a2,nb[i])) or (is_aryl(a2,nb[i])) then inc(c_count);
      if (is_hydroxy(a2,nb[i])) then inc(oh_count);        
      if (is_alkoxy(a2,nb[i])) or (is_aryloxy(a2,nb[i])) then inc(or_count);
      if (is_amino(a2,nb[i])) or (is_subst_amino(a2,nb[i])) then inc(n_count);
      if (nb_el = 'F ') or (nb_el = 'CL') or (nb_el = 'BR') or (nb_el = 'I ') then inc(hal_count);
      if (nb_el = 'O ') then inc(o_count);
    end;  
  if (c_count = 1) then
    begin
      if (atom^[a2].neighbor_count = 3) and
         ((o_count - (oh_count + or_count)) = 1) then    // sulfinic acid & derivs
        begin
          fg[fg_sulfinic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_sulfinic_acid_deriv,a2);   // v0.5
          if (oh_count = 1)  then 
            begin
              fg[fg_sulfinic_acid]        := true;
              if opt_pos then add2fgloc(fg_sulfinic_acid,a2);   // v0.5
            end;
          if (or_count = 1)  then 
            begin
              fg[fg_sulfinic_acid_ester]  := true;
              if opt_pos then add2fgloc(fg_sulfinic_acid_ester,a2);   // v0.5
            end;
          if (hal_count = 1) then 
            begin
              fg[fg_sulfinic_acid_halide] := true;
              if opt_pos then add2fgloc(fg_sulfinic_acid_halide,a2);   // v0.5
            end;
          if (n_count = 1)   then 
            begin
              fg[fg_sulfinic_acid_amide]  := true;
              if opt_pos then add2fgloc(fg_sulfinic_acid_amide,a2);   // v0.5
            end;
        end;
      if (atom^[a2].neighbor_count = 2) and
         ((o_count - (oh_count + or_count)) = 0) then    // sulfenic acid & derivs
        begin
          fg[fg_sulfenic_acid_deriv]  := true;
          if opt_pos then add2fgloc(fg_sulfenic_acid_deriv,a2);   // v0.5
          if (oh_count = 1)  then 
            begin
              fg[fg_sulfenic_acid]        := true;
              if opt_pos then add2fgloc(fg_sulfenic_acid,a2);   // v0.5
            end;
          if (or_count = 1)  then 
            begin
              fg[fg_sulfenic_acid_ester]  := true;
              if opt_pos then add2fgloc(fg_sulfenic_acid_ester,a2);   // v0.5
            end;
          if (hal_count = 1) then 
            begin
              fg[fg_sulfenic_acid_halide] := true;
              if opt_pos then add2fgloc(fg_sulfenic_acid_halide,a2);   // v0.5
            end;
          if (n_count = 1)   then 
            begin
              fg[fg_sulfenic_acid_amide]  := true;
              if opt_pos then add2fgloc(fg_sulfenic_acid_amide,a2);   // v0.5
            end;
        end;
    end;
end;


procedure chk_c_n(a1,a2:integer);
// a1 = C, a2 = N
begin
  // ignore heteroaromatic rings (like furan, thiophene, pyrrol, etc.)
  if (atom^[a2].arom = true) then exit;
  if (is_true_alkyl(a2,a1)) and (is_amino(a1,a2)) then
    begin
      fg[fg_amine]            := true;
      fg[fg_prim_amine]       := true;
      fg[fg_prim_aliph_amine] := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_prim_amine,a2);
          add2fgloc(fg_prim_aliph_amine,a2);
        end;
    end;
  if (is_aryl(a2,a1)) and (is_amino(a1,a2)) then
    begin
      fg[fg_amine]            := true;
      fg[fg_prim_amine]       := true;
      fg[fg_prim_arom_amine]  := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_prim_amine,a2);
          add2fgloc(fg_prim_arom_amine,a2);
        end;

    end;
  if (is_true_alkyl(a2,a1)) and (is_true_alkylamino(a1,a2)) then
    begin
      fg[fg_amine]            := true;
      fg[fg_sec_amine]        := true;
      fg[fg_sec_aliph_amine]  := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_sec_amine,a2);
          add2fgloc(fg_sec_aliph_amine,a2);
        end;
    end;
  if (is_aryl(a2,a1)) and (is_true_alkylamino(a1,a2)) then
    begin
      fg[fg_amine]            := true;
      fg[fg_sec_amine]        := true;
      fg[fg_sec_mixed_amine]  := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_sec_amine,a2);
          add2fgloc(fg_sec_mixed_amine,a2);
        end;
    end;
  if (is_aryl(a2,a1)) and (is_arylamino(a1,a2)) then
    begin
      fg[fg_amine]            := true;
      fg[fg_sec_amine]        := true;
      fg[fg_sec_arom_amine]   := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_sec_amine,a2);
          add2fgloc(fg_sec_arom_amine,a2);
        end;
    end;
  if (is_true_alkyl(a2,a1)) and (is_true_dialkylamino(a1,a2)) then
    begin
      fg[fg_amine]            := true;
      fg[fg_tert_amine]       := true;
      fg[fg_tert_aliph_amine] := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_tert_amine,a2);
          add2fgloc(fg_tert_aliph_amine,a2);
        end;
    end;
  if ((is_true_alkyl(a2,a1)) and (is_diarylamino(a1,a2))) or
     ((is_aryl(a2,a1)) and (is_true_dialkylamino(a1,a2))) then
    begin
      fg[fg_amine]            := true;
      fg[fg_tert_amine]       := true;
      fg[fg_tert_mixed_amine] := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_tert_amine,a2);
          add2fgloc(fg_tert_mixed_amine,a2);
        end;
    end;
  if (is_aryl(a2,a1)) and (is_diarylamino(a1,a2)) then
    begin
      fg[fg_amine]            := true;
      fg[fg_tert_amine]       := true;
      fg[fg_tert_arom_amine]  := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_amine,a2);
          add2fgloc(fg_tert_amine,a2);
          add2fgloc(fg_tert_arom_amine,a2);
        end;
    end;
//  if (is_alkyl(a2,a1) or is_aryl(a2,a1) or    // v0.3k
//      is_alkenyl(a2,a1) or is_alkynyl(a2,a1)) and (is_hydroxylamino(a1,a2) and
//      (is_acyl_gen(a2,a1)=false)) then        // v0.3k
  if ((atom^[a1].element = 'C ') and (is_hydroxylamino(a1,a2)) and (is_acyl_gen(a2,a1)=false) ) then  // v0.5b    
    begin
      fg[fg_hydroxylamine]      := true;        // v0.3k 
      if opt_pos then add2fgloc(fg_hydroxylamine,a2);  // v0.5
    end;
  if (is_alkyl(a2,a1) or is_aryl(a2,a1) or is_acyl(a2,a1) or
      is_alkenyl(a2,a1) or is_alkynyl(a2,a1)) and (is_hydrazino(a1,a2)) then
    begin
      fg[fg_hydrazine]           := true;
      if opt_pos then fgloc_set_hydrazino(a2);
    end;
//  if (is_alkyl(a2,a1) or is_aryl(a2,a1) or    // v0.3k
//      is_alkenyl(a2,a1) or is_alkynyl(a2,a1)) and (is_azido(a1,a2)) then
  if ((atom^[a1].element = 'C ') and (is_azido(a1,a2))) then  // v0.5b    
    begin
      fg[fg_azide]           := true;
      if opt_pos then add2fgloc(fg_azide,a2);   // v0.5
    end;
//  if (is_alkyl(a2,a1) or is_aryl(a2,a1) or    // v0.3k
//      is_alkenyl(a2,a1) or is_alkynyl(a2,a1)) and (is_diazonium(a1,a2)) then
  if ((atom^[a1].element = 'C ') and (is_diazonium(a1,a2))) then  // v0.5b
    begin
      fg[fg_diazonium_salt]           := true;
      if opt_pos then add2fgloc(fg_diazonium_salt,a2);   // v0.5
    end;
//  if (is_alkyl(a2,a1) or is_aryl(a2,a1) or    // v0.3k
//      is_alkenyl(a2,a1) or is_alkynyl(a2,a1)) and (is_nitro(a1,a2)) then
  if ((atom^[a1].element = 'C ') and (is_nitro(a1,a2))) then  // v0.5b
    begin
      fg[fg_nitro_compound]     := true;
      if opt_pos then add2fgloc(fg_nitro_compound,a2);   // v0.5
    end;
  if ((is_alkynyl(a2,a1) and (is_amino(a1,a2) or is_C_monosubst_amino(a1,a2) or 
    is_C_disubst_amino(a1,a2)) and (not is_acylamino(a1,a2)))) then   // v0.4c: fixed parentheses
    begin
      fg[fg_amine]            := true;
      amine_generic           := true;    
      if opt_pos then add2fgloc(fg_amine,a2);   // v0.5
    end;
end;


procedure chk_c_c(a1,a2:integer);
var
  i : integer;
  nb : neighbor_rec;
  oh_count, nhr_count : integer;
  a1oh : boolean;   // v0.5
  b_id : integer;   // v0.5
begin
  // ignore aromatic rings
  if (atom^[a2].arom = true) then exit;
  //check for 1,2-diols and 1,2-aminoalcoholes
  if (atom^[a1].atype = 'C3 ') and (atom^[a2].atype = 'C3 ') then
    begin
      if (hetbond_count(a1) = 1) and (hetbond_count(a2) = 1) then
        begin
          oh_count := 0; nhr_count := 0; a1oh := false;
          fillchar(nb,sizeof(neighbor_rec),0);
          nb := get_neighbors(a1);
          for i := 1 to atom^[a1].neighbor_count do
            begin
              if nb[i] <> a2 then
                begin
                  if (is_hydroxy(a1,nb[i])) then 
                    begin
                      inc(oh_count);        
                      a1oh := true;   // v0.5
                    end;
                  if (is_amino(a1,nb[i])) or (is_alkylamino(a1,nb[i])) 
                      or (is_arylamino(a1,nb[i])) then inc(nhr_count);
                end;
            end;  
          fillchar(nb,sizeof(neighbor_rec),0);
          nb := get_neighbors(a2);
          for i := 1 to atom^[a2].neighbor_count do
            begin
              if nb[i] <> a1 then
                begin
                  if (is_hydroxy(a2,nb[i])) then inc(oh_count);        
                  if (is_amino(a2,nb[i])) or (is_alkylamino(a2,nb[i])) 
                      or (is_arylamino(a2,nb[i])) then inc(nhr_count);
                end;
            end;
          if oh_count = 2 then 
            begin
              fg[fg_1_2_diol] := true;
              if opt_pos then   // v0.5
                begin
                  b_id := get_bond(a1,a2);
                  add2fgloc(fg_1_2_diol,b_id);
                end;
            end;
          if (oh_count = 1) and (nhr_count = 1) then 
            begin
              fg[fg_1_2_aminoalcohol] := true;
              if opt_pos then   // v0.5
                begin
                  if a1oh then add2fgloc(fg_1_2_aminoalcohol,a1) else add2fgloc(fg_1_2_aminoalcohol,a2);
                end;
            end;
        end;
    end;
  // check for alpha-aminoacids and alpha-hydroxyacids
  if (atom^[a1].atype = 'C3 ') and (atom^[a2].atype = 'C2 ') then
    begin
      if (hetbond_count(a1) = 1) and (hetbond_count(a2) = 3) then
        begin
          oh_count := 0; nhr_count := 0;
          fillchar(nb,sizeof(neighbor_rec),0);
          nb := get_neighbors(a1);
          for i := 1 to atom^[a1].neighbor_count do
            begin
              if nb[i] <> a2 then
                begin
                  if (is_hydroxy(a1,nb[i])) then inc(oh_count);        
                  if (is_amino(a1,nb[i])) or (is_alkylamino(a1,nb[i])) 
                      or (is_arylamino(a1,nb[i])) then inc(nhr_count);
                end;
            end;  
          fillchar(nb,sizeof(neighbor_rec),0);
          nb := get_neighbors(a2);
          for i := 1 to atom^[a2].neighbor_count do
            begin
              if nb[i] <> a1 then
                begin
                  if (is_hydroxy(a2,nb[i])) then inc(oh_count);        
                end;
            end;  
          if (oh_count = 2) and (is_oxo_C(a2)) then 
            begin
              fg[fg_alpha_hydroxyacid] := true;
              if opt_pos then add2fgloc(fg_alpha_hydroxyacid,a2);   // v0.5
            end;
          if (oh_count = 1) and (nhr_count = 1) and (is_oxo_C(a2)) then 
            begin
              fg[fg_alpha_aminoacid] := true;
              if opt_pos then add2fgloc(fg_alpha_aminoacid,a2);   // v0.5
            end;
        end;
    end;    
end;


procedure chk_x_y_single(a_view,a_ref:integer);
var
  b_id : integer;  // v0.5
begin
  if (atom^[a_view].atype = 'O3 ') and (atom^[a_ref].atype = 'O3 ') then
    begin
      if (is_hydroxy(a_ref,a_view)) or
         (is_hydroxy(a_view,a_ref)) then 
        begin
          fg[fg_hydroperoxide] := true;
          if opt_pos then   // v0.5
            begin
              if is_hydroxy(a_view,a_ref) then add2fgloc(fg_hydroperoxide,a_view) else
                add2fgloc(fg_hydroperoxide,a_ref);
            end;
        end;
      if ((is_alkoxy(a_ref,a_view)) or
          (is_aryloxy(a_ref,a_view)) or 
          (is_siloxy(a_ref,a_view))) and
         ((is_alkoxy(a_view,a_ref)) or 
          (is_aryloxy(a_view,a_ref)) or
          (is_siloxy(a_view,a_ref))) then 
        begin
          fg[fg_peroxide] := true;
          if opt_pos then
            begin
              b_id := get_bond(a_view,a_ref);
              add2fgloc(fg_peroxide,b_id);
            end;
        end;
    end;  // still missing: peracid
  if (atom^[a_view].atype = 'S3 ') and (atom^[a_ref].atype = 'S3 ') then
    begin
      if (atom^[a_view].neighbor_count = 2) and (atom^[a_ref].neighbor_count = 2) then
        begin
          fg[fg_disulfide] := true;
          if opt_pos then   // v0.5
            begin
              b_id := get_bond(a_view,a_ref);
              add2fgloc(fg_disulfide,b_id);
            end;
        end;
    end;
  if (atom^[a_view].element = 'N ') and (atom^[a_ref].element = 'N ') and
     (hetbond_count(a_view) = 1) and (hetbond_count(a_ref) = 1) then
    begin
      //if ((is_amino(a_ref,a_view)) or 
      //    (is_subst_amino(a_ref,a_view)) or
      //    (is_acylamino(a_ref,a_view))) and
      //   ((is_amino(a_view,a_ref)) or 
      //    (is_subst_amino(a_view,a_ref)) or
      //    (is_acylamino(a_ref,a_view))) then 
      if (bond^[get_bond(a_view,a_ref)].arom = false) then 
        begin
          fg[fg_hydrazine] := true;
          if opt_pos then   // v0.5
            begin
              b_id := get_bond(a_view,a_ref);
              add2fgloc(fg_hydrazine,b_id);
            end;

        end;
    end; 
  if (atom^[a_view].element = 'N ') and (atom^[a_ref].atype = 'O3 ') then
    begin  // bond is in "opposite" direction
      if ((is_alkoxy(a_view,a_ref)) or (is_aryloxy(a_view,a_ref))) and
         (is_nitro(a_ref,a_view)) then 
        begin
          fg[fg_nitrate] := true;
          if opt_pos then add2fgloc(fg_nitrate,a_view);   // v0.5
        end;
      if ((is_nitro(a_ref,a_view)=false) and (atom^[a_view].arom=false)) and
         ((is_amino(a_ref,a_view)) or (is_subst_amino(a_ref,a_view))) and
         (is_acylamino(a_ref,a_view)=false) then
        begin
          fg[fg_hydroxylamine] := true;    // new in v0.3c
          if opt_pos then add2fgloc(fg_hydroxylamine,a_view);   // v0.5
        end;
    end;
  if (atom^[a_view].element = 'S ') and (atom^[a_ref].element = 'O ') then    
    chk_sulfoxide(a_view,a_ref);
end;


procedure chk_single(a1,a2:integer);
var
  a1_el, a2_el : str2;
begin
  a1_el := atom^[a1].element;
  a2_el := atom^[a2].element;
  if (a1_el = 'C ') and 
     ((a2_el = 'F ') or (a2_el = 'CL') or (a2_el = 'BR') or (a2_el = 'I ')) then chk_c_hal(a1,a2);
  if (a1_el = 'C ') and (a2_el = 'O ') then chk_c_o(a1,a2);   
  if (a1_el = 'C ') and (a2_el = 'S ') then chk_c_s(a1,a2);   
  if (a1_el = 'C ') and (a2_el = 'N ') then chk_c_n(a1,a2);       
  if (a1_el = 'C ') and (atom^[a2].metal and (is_cyano_c(a1) = false)) then
    begin
      fg[fg_organometallic] := true;
      if opt_pos then add2fgloc(fg_organometallic,a2);   // v0.5
      if (a2_el = 'LI') then 
        begin
          fg[fg_organolithium] := true;
          if opt_pos then add2fgloc(fg_organolithium,a2);   // v0.5
        end;
      if (a2_el = 'MG') then 
        begin
           fg[fg_organomagnesium] := true;
           if opt_pos then add2fgloc(fg_organomagnesium,a2);   // v0.5
         end;
    end; 
  if (a1_el = 'C ') and (a2_el = 'C ') then chk_c_c(a1,a2);       
  if (a1_el <> 'C ') and (a2_el <> 'C ') then chk_x_y_single(a1,a2);
end;


procedure chk_carbonyl_deriv_sp3(a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  oh_count, or_count, n_count, sh_count, sr_count : integer;
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_ref);
  oh_count := 0; or_count := 0; n_count := 0; sh_count := 0; sr_count := 0;
  for i := 1 to atom^[a_ref].neighbor_count do
    begin
      if (is_hydroxy(a_ref,nb[i])) then inc(oh_count);
      if (is_alkoxy(a_ref,nb[i])) or (is_aryloxy(a_ref,nb[i])) or
         (is_alkenyloxy(a_ref,nb[i])) or (is_alkynyloxy(a_ref,nb[i])) then inc(or_count);
      if (is_sulfanyl(a_ref,nb[i])) then inc(sh_count);
      if (is_alkylsulfanyl(a_ref,nb[i])) or (is_arylsulfanyl(a_ref,nb[i])) or
         (is_alkenylsulfanyl(a_ref,nb[i])) or (is_alkynylsulfanyl(a_ref,nb[i])) then inc(sr_count);
      if (atom^[(nb[i])].atype = 'N3 ') or (atom^[(nb[i])].atype = 'NAM') then inc(n_count);
    end;
  if (oh_count = 2) then 
    begin
      fg[fg_carbonyl_hydrate] := true;
      if opt_pos then add2fgloc(fg_carbonyl_hydrate,a_ref);  // v0.5
    end;
  if (oh_count = 1) and (or_count = 1) then 
    begin
      fg[fg_hemiacetal] := true;
      if opt_pos then add2fgloc(fg_hemiacetal,a_ref);  // v0.5
    end;
  if (or_count = 2) then 
    begin
      fg[fg_acetal] := true;
      if opt_pos then add2fgloc(fg_acetal,a_ref);  // v0.5
    end;
  if ((oh_count = 1) or (or_count = 1)) and (n_count = 1) then 
    begin
      fg[fg_hemiaminal] := true;  
      if opt_pos then add2fgloc(fg_hemiaminal,a_ref);  // v0.5
    end;
  if (n_count = 2) then 
    begin
      fg[fg_aminal] := true;  
      if opt_pos then add2fgloc(fg_aminal,a_ref);  // v0.5
    end;
  if ((sh_count = 1) or (sr_count = 1)) and (n_count = 1) then 
    begin
      fg[fg_thiohemiaminal] := true;  
      if opt_pos then add2fgloc(fg_thiohemiaminal,a_ref);  // v0.5
    end;
  if (sr_count = 2) or ((or_count = 1) and (sr_count = 1)) then 
    begin
      fg[fg_thioacetal] := true;  
      if opt_pos then add2fgloc(fg_thioacetal,a_ref);  // v0.5
    end;
end;


procedure chk_carboxyl_deriv_sp3(a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  or_count, oh_count, n_count : integer;  // oh_count new in v0.3c
  electroneg_count : integer;             // new in v0.3j
  hal_count        : integer;
  nb_el : str2;
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_ref);
  or_count := 0; oh_count := 0; n_count := 0;
  electroneg_count := 0;
  hal_count        := 0;
  for i := 1 to atom^[a_ref].neighbor_count do
    begin
      nb_el := atom^[(nb[i])].element;  // v0.3j
      if (is_electroneg(nb_el)) then inc(electroneg_count); 
      if ((nb_el='F ') or (nb_el='CL') or (nb_el='BR') or (nb_el='I ')) then inc(hal_count);
      if (is_alkoxy(a_ref,nb[i])) or
         (is_aryloxy(a_ref,nb[i])) or
         (is_siloxy(a_ref,nb[i])) then inc(or_count);
      if (is_hydroxy(a_ref,nb[i])) then inc(oh_count);  // new in v0.3c   
      if (atom^[(nb[i])].atype = 'N3 ') or (atom^[(nb[i])].atype = 'NAM') then inc(n_count);
    end;
  //if (or_count + n_count > 1) then fg[fg_orthocarboxylic_acid_deriv] := true;  // until v0.3i
  if (electroneg_count = 3) and (hal_count < 3) then 
    begin
      fg[fg_orthocarboxylic_acid_deriv] := true;  // v0.3j
      if opt_pos then add2fgloc(fg_orthocarboxylic_acid_deriv,a_ref);   // v0.5
    end;
  if (or_count = 3) then 
    begin
      fg[fg_carboxylic_acid_orthoester] := true;
      if opt_pos then add2fgloc(fg_carboxylic_acid_orthoester,a_ref);   // v0.5
    end;
  if (or_count = 2) and (n_count = 1) then 
    begin
      fg[fg_carboxylic_acid_amide_acetal] := true;
      if opt_pos then add2fgloc(fg_carboxylic_acid_amide_acetal,a_ref);   // v0.5
    end;
  if (oh_count > 0) and (oh_count + or_count + n_count = 3) then 
    begin
      fg[fg_orthocarboxylic_acid_deriv] := true;  // new in v0.3c
      if opt_pos then add2fgloc(fg_orthocarboxylic_acid_deriv,a_ref);   // v0.5
    end;
end;


procedure chk_anhydride(a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  acyl_count : integer;
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_ref);
  acyl_count := 0;
  for i := 1 to atom^[a_ref].neighbor_count do
    begin
      if is_acyl_gen(a_ref,nb[i]) then inc(acyl_count); // v0.4b replaced is_acyl() by is_acyl_gen()
    end;
  if (acyl_count = 2) and (atom^[a_ref].atype = 'O3 ') then
    begin
      fg[fg_carboxylic_acid_deriv]     := true;
      fg[fg_carboxylic_acid_anhydride] := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_carboxylic_acid_anhydride,a_ref);
          for i := 1 to atom^[a_ref].neighbor_count do
            begin
              if is_acyl(a_ref,nb[i]) then add2fgloc(fg_carboxylic_acid_deriv,nb[i]);  // use is_acyl, not is_acyl_gen here!
            end;
        end;
    end;
end;


procedure chk_imide(a_ref:integer);
var
  i : integer;
  nb : neighbor_rec;
  acyl_count : integer;
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_ref);
  acyl_count := 0;
  for i := 1 to atom^[a_ref].neighbor_count do
    begin
      if (is_acyl_gen(a_ref,nb[i])) or (is_carbamoyl(a_ref,nb[i])) then inc(acyl_count);  // v0.3j
    end;
  if (acyl_count >= 2) and (atom^[a_ref].element = 'N ') then  // v0.3j: accept also N-acyl-imides
    begin
      fg[fg_carboxylic_acid_deriv]     := true;
      fg[fg_carboxylic_acid_imide] := true;
      if (atom^[a_ref].neighbor_count = 2) then 
        fg[fg_carboxylic_acid_unsubst_imide] := true;
      if (atom^[a_ref].neighbor_count = 3) then 
        fg[fg_carboxylic_acid_subst_imide] := true;
      if opt_pos then   // v0.5
        begin
          add2fgloc(fg_carboxylic_acid_imide,a_ref);
          if (atom^[a_ref].neighbor_count = 2) then add2fgloc(fg_carboxylic_acid_unsubst_imide,a_ref);
          if (atom^[a_ref].neighbor_count = 3) then  add2fgloc(fg_carboxylic_acid_subst_imide,a_ref);
          for i := 1 to atom^[a_ref].neighbor_count do
            begin
              if is_acyl(a_ref,nb[i]) then add2fgloc(fg_carboxylic_acid_deriv,nb[i]);  // use is_acyl, not is_acyl_gen here!
            end;
        end;
    end;
end;


procedure chk_12diphenol(a_view,a_ref:integer);
// a_view and a_ref are adjacent ring atoms
var
  i : integer;
  nb : neighbor_rec;
  oh_count : integer;
  b_id : integer;  // v0.5
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_view);
  oh_count := 0;
  for i := 1 to atom^[a_view].neighbor_count do
    begin
      if bond^[(get_bond(a_view,nb[i]))].btype = 'S' then
        begin
          if (is_hydroxy(a_view,nb[i])) then inc(oh_count);
        end;
    end;
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_ref);
  for i := 1 to atom^[a_ref].neighbor_count do
    begin
      if bond^[(get_bond(a_ref,nb[i]))].btype = 'S' then
        begin
          if (is_hydroxy(a_ref,nb[i])) then inc(oh_count);
        end;
    end;
  if (oh_count = 2) then 
    begin
      fg[fg_1_2_diphenol] := true;
      if opt_pos then   // v0.5
        begin
          b_id := get_bond(a_view,a_ref);
          add2fgloc(fg_1_2_diphenol,b_id);
        end;
    end;
end;


procedure chk_arom_fg(a1,a2:integer);
begin
  if (hetbond_count(a1) = 1) and (hetbond_count(a2) = 1) then chk_12diphenol(a1,a2);
end;

function is_arene(r_id:integer):boolean;
var
  i,j  : integer;
  r    : boolean;
  testring : ringpath_type;
  ring_size : integer;
  a_prev, a_ref : integer;
begin
  r := false;
  if (r_id < 1) or (r_id > n_rings) then begin is_arene := false; exit; end;
  fillchar(testring,sizeof(ringpath_type),0);
  ring_size := ringprop^[r_id].size;  // v0.3j
  //for j := 1 to max_ringsize do if ring^[r_id,j] > 0 then testring[j] := ring^[r_id,j];
  for j := 1 to ring_size do testring[j] := ring^[r_id,j];  // v0.3j
  //ring_size := path_length(testring);
  if (ring_size > 2) then
    begin
      r := true;
      a_prev := testring[ring_size];
      for i := 1 to ring_size do
        begin
          a_ref := testring[i];
          if (bond^[get_bond(a_prev,a_ref)].arom = false) then r := false;
          a_prev := a_ref;
        end;
    end;
  is_arene := r;
end;


function is_heterocycle(r_id:integer):boolean;
var
  i,j  : integer;
  r    : boolean;
  testring : ringpath_type;
  ring_size : integer;
  a_ref : integer;
begin
  r := false;
  if (r_id < 1) or (r_id > n_rings) then begin is_heterocycle := false; exit; end;
  fillchar(testring,sizeof(ringpath_type),0);
  ring_size := ringprop^[r_id].size;  // v0.3j
  //for j := 1 to max_ringsize do if ring^[r_id,j] > 0 then testring[j] := ring^[r_id,j];
  for j := 1 to ring_size do testring[j] := ring^[r_id,j];  // v0.3j
  //ring_size := path_length(testring);
  if (ring_size > 2) then
    begin
      for i := 1 to ring_size do
        begin
          a_ref := testring[i];
          if (atom^[a_ref].element <> 'C ') then r := true;
        end;
    end;
  is_heterocycle := r;
end;


procedure chk_oxo_thioxo_imino_hetarene(r_id:integer);
var
  i,j  : integer;
  testring : ringpath_type;
  ring_size : integer;
  a_ref : integer;
begin
  if (r_id < 1) or (r_id > n_rings) then exit;
  fillchar(testring,sizeof(ringpath_type),0);
  ring_size := ringprop^[r_id].size;  // v0.3j
  //for j := 1 to max_ringsize do if ring^[r_id,j] > 0 then testring[j] := ring^[r_id,j];
  for j := 1 to ring_size do testring[j] := ring^[r_id,j];  // v0.3j
  //ring_size := path_length(testring);
  //if (is_arene(r_id)) and (odd(ring_size) = false) then
  if is_arene(r_id) then    // v0.3j
    begin
      for i := 1 to ring_size do
        begin
          a_ref := testring[i];
          if (is_oxo_c(a_ref)) then 
            begin
              fg[fg_oxohetarene]  := true;
              if opt_pos then add2fgloc(fg_oxohetarene,a_ref);   // v0.5
            end;
          if (is_thioxo_c(a_ref)) then 
            begin
              fg[fg_thioxohetarene]  := true;
              if opt_pos then add2fgloc(fg_thioxohetarene,a_ref);   // v0.5
            end;
          if (is_true_exocyclic_imino_c(a_ref,r_id)) then   // v0.3j
            begin
              fg[fg_iminohetarene] := true;
              if opt_pos then add2fgloc(fg_iminohetarene,a_ref);   // v0.5
            end;
        end;
    end;
end;


procedure chk_ion(a_ref:integer);
var
  i      : integer;
  nb     : neighbor_rec;
  charge : integer;
begin
  fillchar(nb,sizeof(neighbor_rec),0);
  nb := get_neighbors(a_ref);
  charge := atom^[a_ref].formal_charge;
  if charge <> 0 then
    begin  // check if charge is neutralized by an adjacent opposite charge
      for i := 1 to atom^[a_ref].neighbor_count do
        begin
          charge := charge + atom^[(nb[i])].formal_charge;
        end;
      if (charge > 0) then 
        begin
          fg[fg_cation] := true;
          if opt_pos then add2fgloc(fg_cation,a_ref);  // v0.5
        end;
      if (charge < 0) then 
        begin
          fg[fg_anion]  := true;
          if opt_pos then add2fgloc(fg_anion,a_ref);  // v0.5
        end;
    end;
end;


procedure chk_functionalgroups;
var
  i : integer;
  a1,a2 : integer;
  bt : char;
  pos_chg, neg_chg : integer;
begin
  if (n_atoms < 1) or (n_bonds < 1) then exit;
  if opt_pos then  // v0.5
    begin  // allocate memory for the location table
      if (fgloc = nil) then 
        try
          getmem(fgloc,sizeof(fgloctype));
        except
          on e:Eoutofmemory do
            begin
              writeln('Not enough memory for -p option');
              opt_pos := false;
            end;
        end;
      fillchar(fgloc^,sizeof(fgloctype),0);
    end;
  pos_chg := 0; neg_chg := 0;
  for i := 1 to n_atoms do  // a few groups are best discovered in the atom list
    begin
      if (atom^[i].atype = 'SO2') then chk_so2_deriv(i);
      //if (atom^[i].atype = 'SO ') then fg[fg_sulfoxide] := true;  // do another check in the bond list!!
      if (atom^[i].element = 'P ') then chk_p_deriv(i);
      if (atom^[i].element = 'B ') then chk_b_deriv(i);
      if (atom^[i].atype = 'N3+') or (atom^[i].formal_charge > 0) then chk_ammon(i);
      if (atom^[i].atype = 'C3 ') and (hetbond_count(i) = 2) then chk_carbonyl_deriv_sp3(i);
      if (atom^[i].atype = 'C3 ') and (hetbond_count(i) = 3) then chk_carboxyl_deriv_sp3(i);
      if (atom^[i].atype = 'O3 ') and (atom^[i].neighbor_count = 2) then chk_anhydride(i);
      if ((atom^[i].atype = 'N3 ') or (atom^[i].atype = 'NAM')) and
          (atom^[i].neighbor_count >= 2) then chk_imide(i);
      if (atom^[i].formal_charge > 0) then pos_chg := pos_chg + atom^[i].formal_charge;
      if (atom^[i].formal_charge < 0) then neg_chg := neg_chg + atom^[i].formal_charge;
      chk_ion(i);
    end;
  for i := 1 to n_bonds do  // most groups are best discovered in the bond list
    begin
      a1 := bond^[i].a1;
      a2 := bond^[i].a2;
      bt := bond^[i].btype;
      if (atom^[a1].heavy) and (atom^[a2].heavy) then
        begin
          orient_bond(a1,a2);
          if (bt = 'T') then chk_triple(a1,a2);
          if (bt = 'D') then chk_double(a1,a2);
          if (bt = 'S') then chk_single(a1,a2);
          if (bond^[i].arom) then chk_arom_fg(a1,a2);
        end;
    end;
  if (n_rings > 0) then
    begin
      for i := 1 to n_rings do
        begin
          chk_oxo_thioxo_imino_hetarene(i);
          if (is_arene(i)) then 
            begin
              fg[fg_aromatic] := true;
              if (opt_pos and (ringprop^[i].envelope = false)) then add2fgloc(fg_aromatic,i);    // v0.5
            end;
          if (is_heterocycle(i)) then 
            begin
              fg[fg_heterocycle] := true;
              if (opt_pos and (ringprop^[i].envelope = false)) then add2fgloc(fg_heterocycle,i);    // v0.5
            end;
        end;
    end;
  if (pos_chg + neg_chg) > 0 then fg[fg_cation] := true;
  if (pos_chg + neg_chg) < 0 then fg[fg_anion]  := true;
end;

(*
procedure write_fg_text;   // old version
begin
  if fg[fg_cation]                         then writeln('cation');
  if fg[fg_anion]                          then writeln('anion');
//  if fg[fg_carbonyl]                       then writeln('carbonyl compound');
  if fg[fg_aldehyde]                       then writeln('aldehyde');
  if fg[fg_ketone]                         then writeln('ketone');
//  if fg[fg_thiocarbonyl]                   then writeln('thiocarbonyl compound');
  if fg[fg_thioaldehyde]                   then writeln('thioaldehyde');
  if fg[fg_thioketone]                     then writeln('thioketone');
  if fg[fg_imine]                          then writeln('imine');
  if fg[fg_hydrazone]                      then writeln('hydrazone');
  if fg[fg_semicarbazone]                  then writeln('semicarbazone');
  if fg[fg_thiosemicarbazone]              then writeln('thiosemicarbazone');
  if fg[fg_oxime]                          then writeln('oxime');
  if fg[fg_oxime_ether]                    then writeln('oxime ether');
  if fg[fg_ketene]                         then writeln('ketene');
  if fg[fg_ketene_acetal_deriv]            then writeln('ketene acetal or derivative');
  if fg[fg_carbonyl_hydrate]               then writeln('carbonyl hydrate');
  if fg[fg_hemiacetal]                     then writeln('hemiacetal');
  if fg[fg_acetal]                         then writeln('acetal');
  if fg[fg_hemiaminal]                     then writeln('hemiaminal');
  if fg[fg_aminal]                         then writeln('aminal');
  if fg[fg_thiohemiaminal]                 then writeln('hemithioaminal');
  if fg[fg_thioacetal]                     then writeln('thioacetal');
  if fg[fg_enamine]                        then writeln('enamine');
  if fg[fg_enol]                           then writeln('enol');
  if fg[fg_enolether]                      then writeln('enol ether');
  if fg[fg_hydroxy] and hydroxy_generic    then writeln('hydroxy compound');
//  if fg[fg_alcohol]                        then writeln('alcohol');
  if fg[fg_prim_alcohol]                   then writeln('primary alcohol');
  if fg[fg_sec_alcohol]                    then writeln('secondary alcohol');
  if fg[fg_tert_alcohol]                   then writeln('tertiary alcohol');
  if fg[fg_1_2_diol]                       then writeln('1,2-diol');
  if fg[fg_1_2_aminoalcohol]               then writeln('1,2-aminoalcohol');
  if fg[fg_phenol]                         then writeln('phenol or hydroxyhetarene');
  if fg[fg_1_2_diphenol]                   then writeln('1,2-diphenol');
  if fg[fg_enediol]                        then writeln('enediol');
  if fg[fg_ether] and ether_generic        then writeln('ether');
  if fg[fg_dialkylether]                   then writeln('dialkyl ether');
  if fg[fg_alkylarylether]                 then writeln('alkyl aryl ether ');
  if fg[fg_diarylether]                    then writeln('diaryl ether');
  if fg[fg_thioether]                      then writeln('thioether');
  if fg[fg_disulfide]                      then writeln('disulfide');
  if fg[fg_peroxide]                       then writeln('peroxide');
  if fg[fg_hydroperoxide]                  then writeln('hydroperoxide ');
  if fg[fg_hydrazine]                      then writeln('hydrazine derivative');
  if fg[fg_hydroxylamine]                  then writeln('hydroxylamine');
  if fg[fg_amine] and amine_generic        then writeln('amine');
  if fg[fg_prim_amine]                     then writeln('primary amine');
  if fg[fg_prim_aliph_amine]               then writeln('primary aliphatic amine (alkylamine)');
  if fg[fg_prim_arom_amine]                then writeln('primary aromatic amine');
  if fg[fg_sec_amine]                      then writeln('secondary amine');
  if fg[fg_sec_aliph_amine]                then writeln('secondary aliphatic amine (dialkylamine)');
  if fg[fg_sec_mixed_amine]                then writeln('secondary aliphatic/aromatic amine (alkylarylamine)');
  if fg[fg_sec_arom_amine]                 then writeln('secondary aromatic amine (diarylamine)');
  if fg[fg_tert_amine]                     then writeln('tertiary amine');
  if fg[fg_tert_aliph_amine]               then writeln('tertiary aliphatic amine (trialkylamine)');
  if fg[fg_tert_mixed_amine]               then writeln('tertiary aliphatic/aromatic amine (alkylarylamine)');
  if fg[fg_tert_arom_amine]                then writeln('tertiary aromatic amine (triarylamine)');
  if fg[fg_quart_ammonium]                 then writeln('quaternary ammonium salt');
  if fg[fg_n_oxide]                        then writeln('N-oxide');
  // new in v0.2f
  if fg[fg_halogen_deriv]                  then 
    begin
      if (not fg[fg_alkyl_halide]) and (not fg[fg_aryl_halide]) and (not fg[fg_acyl_halide]) then
      writeln('halogen derivative');
    end;
//  if fg[fg_alkyl_halide]                   then writeln('alkyl halide');
  if fg[fg_alkyl_fluoride]                 then writeln('alkyl fluoride');
  if fg[fg_alkyl_chloride]                 then writeln('alkyl chloride');
  if fg[fg_alkyl_bromide]                  then writeln('alkyl bromide');
  if fg[fg_alkyl_iodide]                   then writeln('alkyl iodide');
//  if fg[fg_aryl_halide]                    then writeln('aryl halide');
  if fg[fg_aryl_fluoride]                  then writeln('aryl fluoride');
  if fg[fg_aryl_chloride]                  then writeln('aryl chloride');
  if fg[fg_aryl_bromide]                   then writeln('aryl bromide');
  if fg[fg_aryl_iodide]                    then writeln('aryl iodide');
  if fg[fg_organometallic]                 then writeln('organometallic compound');
  if fg[fg_organolithium]                  then writeln('organolithium compound');
  if fg[fg_organomagnesium]                then writeln('organomagnesium compound');
//  if fg[fg_carboxylic_acid_deriv]          then writeln('carboxylic acid derivative');
  if fg[fg_carboxylic_acid]                then writeln('carboxylic acid');
  if fg[fg_carboxylic_acid_salt]           then writeln('carboxylic acid salt');
  if fg[fg_carboxylic_acid_ester]          then writeln('carboxylic acid ester');
  if fg[fg_lactone]                        then writeln('lactone');
//  if fg[fg_carboxylic_acid_amide]          then writeln('carboxylic acid amide');
  if fg[fg_carboxylic_acid_prim_amide]     then writeln('primary carboxylic acid amide');
  if fg[fg_carboxylic_acid_sec_amide]      then writeln('secondary carboxylic acid amide');
  if fg[fg_carboxylic_acid_tert_amide]     then writeln('tertiary carboxylic acid amide');
  if fg[fg_lactam]                         then writeln('lactam');
  if fg[fg_carboxylic_acid_hydrazide]      then writeln('carboxylic acid hydrazide');
  if fg[fg_carboxylic_acid_azide]          then writeln('carboxylic acid azide');
  if fg[fg_hydroxamic_acid]                then writeln('hydroxamic acid');
  if fg[fg_carboxylic_acid_amidine]        then writeln('carboxylic acid amidine');
  if fg[fg_carboxylic_acid_amidrazone]     then writeln('carboxylic acid amidrazone');
  if fg[fg_nitrile]                        then writeln('carbonitrile');
//  if fg[fg_acyl_halide]                    then writeln('acyl halide');
  if fg[fg_acyl_fluoride]                  then writeln('acyl fluoride');
  if fg[fg_acyl_chloride]                  then writeln('acyl chloride');
  if fg[fg_acyl_bromide]                   then writeln('acyl bromide');
  if fg[fg_acyl_iodide]                    then writeln('acyl iodide');
  if fg[fg_acyl_cyanide]                   then writeln('acyl cyanide');
  if fg[fg_imido_ester]                    then writeln('imido ester');
  if fg[fg_imidoyl_halide]                 then writeln('imidoyl halide');
//  if fg[fg_thiocarboxylic_acid_deriv]      then writeln('thiocarboxylic acid derivative');
  if fg[fg_thiocarboxylic_acid]            then writeln('thiocarboxylic acid');
  if fg[fg_thiocarboxylic_acid_ester]      then writeln('thiocarboxylic acid ester');
  if fg[fg_thiolactone]                    then writeln('thiolactone');
  if fg[fg_thiocarboxylic_acid_amide]      then writeln('thiocarboxylic acid amide');
  if fg[fg_thiolactam]                     then writeln('thiolactam');
  if fg[fg_imido_thioester]                then writeln('imidothioester');
  if fg[fg_oxohetarene]                    then writeln('oxo(het)arene');
  if fg[fg_thioxohetarene]                 then writeln('thioxo(het)arene');
  if fg[fg_iminohetarene]                  then writeln('imino(het)arene');
  if fg[fg_orthocarboxylic_acid_deriv]     then writeln('orthocarboxylic acid derivative');
  if fg[fg_carboxylic_acid_orthoester]     then writeln('orthoester');
  if fg[fg_carboxylic_acid_amide_acetal]   then writeln('amide acetal');
  if fg[fg_carboxylic_acid_anhydride]      then writeln('carboxylic acid anhydride');
//  if fg[fg_carboxylic_acid_imide]          then writeln('carboxylic acid imide');
  if fg[fg_carboxylic_acid_unsubst_imide]  then writeln('carboxylic acid imide, N-unsubstituted');
  if fg[fg_carboxylic_acid_subst_imide]    then writeln('carboxylic acid imide, N-substituted');
  if fg[fg_co2_deriv]                      then writeln('CO2 derivative (general)');
  if fg[fg_carbonic_acid_deriv] and not    // changed in v0.3c
    (fg[fg_carbonic_acid_monoester] or 
     fg[fg_carbonic_acid_diester] or
     fg[fg_carbonic_acid_ester_halide])    then writeln('carbonic acid derivative');
  if fg[fg_carbonic_acid_monoester]        then writeln('carbonic acid monoester');
  if fg[fg_carbonic_acid_diester]          then writeln('carbonic acid diester');
  if fg[fg_carbonic_acid_ester_halide]     then writeln('carbonic acid ester halide (alkyl/aryl haloformate)');
  if fg[fg_thiocarbonic_acid_deriv]        then writeln('thiocarbonic acid derivative');
  if fg[fg_thiocarbonic_acid_monoester]    then writeln('thiocarbonic acid monoester');
  if fg[fg_thiocarbonic_acid_diester]      then writeln('thiocarbonic acid diester');
  if fg[fg_thiocarbonic_acid_ester_halide] then writeln('thiocarbonic acid ester halide (alkyl/aryl halothioformate)');
  if fg[fg_carbamic_acid_deriv] and not    // changed in v0.3c
    (fg[fg_carbamic_acid] or
     fg[fg_carbamic_acid_ester] or
     fg[fg_carbamic_acid_halide])          then writeln('carbamic acid derivative');
  if fg[fg_carbamic_acid]                  then writeln('carbamic acid');
  if fg[fg_carbamic_acid_ester]            then writeln('carbamic acid ester (urethane)');
  if fg[fg_carbamic_acid_halide]           then writeln('carbamic acid halide (haloformic acid amide)');
  if fg[fg_thiocarbamic_acid_deriv] and not  // changed in v0.3c
     (fg[fg_thiocarbamic_acid] or
      fg[fg_thiocarbamic_acid_ester] or
      fg[fg_thiocarbamic_acid_halide])     then writeln('thiocarbamic acid derivative');
  if fg[fg_thiocarbamic_acid]              then writeln('thiocarbamic acid');
  if fg[fg_thiocarbamic_acid_ester]        then writeln('thiocarbamic acid ester');
  if fg[fg_thiocarbamic_acid_halide]       then writeln('thiocarbamic acid halide (halothioformic acid amide)');
  if fg[fg_urea]                           then writeln('urea');
  if fg[fg_isourea]                        then writeln('isourea');
  if fg[fg_thiourea]                       then writeln('thiourea');
  if fg[fg_isothiourea]                    then writeln('isothiourea');
  if fg[fg_guanidine]                      then writeln('guanidine');
  if fg[fg_semicarbazide]                  then writeln('semicarbazide');
  if fg[fg_thiosemicarbazide]              then writeln('thiosemicarbazide');
  if fg[fg_azide]                          then writeln('azide');
  if fg[fg_azo_compound]                   then writeln('azo compound');
  if fg[fg_diazonium_salt]                 then writeln('diazonium salt');
  if fg[fg_isonitrile]                     then writeln('isonitrile');
  if fg[fg_cyanate]                        then writeln('cyanate');
  if fg[fg_isocyanate]                     then writeln('isocyanate');
  if fg[fg_thiocyanate]                    then writeln('thiocyanate');
  if fg[fg_isothiocyanate]                 then writeln('isothiocyanate');
  if fg[fg_carbodiimide]                   then writeln('carbodiimide');
  if fg[fg_nitroso_compound]               then writeln('nitroso compound');
  if fg[fg_nitro_compound]                 then writeln('nitro compound');
  if fg[fg_nitrite]                        then writeln('nitrite');
  if fg[fg_nitrate]                        then writeln('nitrate');
//  if fg[fg_sulfuric_acid_deriv]            then writeln('sulfuric acid derivative');
  if fg[fg_sulfuric_acid]                  then writeln('sulfuric acid');
  if fg[fg_sulfuric_acid_monoester]        then writeln('sulfuric acid monoester');
  if fg[fg_sulfuric_acid_diester]          then writeln('sulfuric acid diester');
  if fg[fg_sulfuric_acid_amide_ester]      then writeln('sulfuric acid amide ester');
  if fg[fg_sulfuric_acid_amide]            then writeln('sulfuric acid amide');
  if fg[fg_sulfuric_acid_diamide]          then writeln('sulfuric acid diamide');
  if fg[fg_sulfuryl_halide]                then writeln('sulfuryl halide');
//  if fg[fg_sulfonic_acid_deriv]            then writeln('sulfonic acid derivative ');
  if fg[fg_sulfonic_acid]                  then writeln('sulfonic acid');
  if fg[fg_sulfonic_acid_ester]            then writeln('sulfonic acid ester');
  if fg[fg_sulfonamide]                    then writeln('sulfonamide');
  if fg[fg_sulfonyl_halide]                then writeln('sulfonyl halide');
  if fg[fg_sulfone]                        then writeln('sulfone');
  if fg[fg_sulfoxide]                      then writeln('sulfoxide');
//  if fg[fg_sulfinic_acid_deriv]            then writeln('sulfinic acid derivative');
  if fg[fg_sulfinic_acid]                  then writeln('sulfinic acid');
  if fg[fg_sulfinic_acid_ester]            then writeln('sulfinic acid ester');
  if fg[fg_sulfinic_acid_halide]           then writeln('sulfinic acid halide');
  if fg[fg_sulfinic_acid_amide]            then writeln('sulfinic acid amide');
//  if fg[fg_sulfenic_acid_deriv]            then writeln('sulfenic acid derivative');
  if fg[fg_sulfenic_acid]                  then writeln('sulfenic acid');
  if fg[fg_sulfenic_acid_ester]            then writeln('sulfenic acid ester');
  if fg[fg_sulfenic_acid_halide]           then writeln('sulfenic acid halide');
  if fg[fg_sulfenic_acid_amide]            then writeln('sulfenic acid amide');
  if fg[fg_thiol]                          then writeln('thiol (sulfanyl compound)');
  if fg[fg_alkylthiol]                     then writeln('alkylthiol');
  if fg[fg_arylthiol]                      then writeln('arylthiol');
//  if fg[fg_phosphoric_acid_deriv]          then writeln('phosphoric acid derivative');
  if fg[fg_phosphoric_acid]                then writeln('phosphoric acid');
  if fg[fg_phosphoric_acid_ester]          then writeln('phosphoric acid ester');
  if fg[fg_phosphoric_acid_halide]         then writeln('phosphoric acid halide');
  if fg[fg_phosphoric_acid_amide]          then writeln('phosphoric acid amide');
//  if fg[fg_thiophosphoric_acid_deriv]      then writeln('thiophosphoric acid derivative');
  if fg[fg_thiophosphoric_acid]            then writeln('thiophosphoric acid');
  if fg[fg_thiophosphoric_acid_ester]      then writeln('thiophosphoric acid ester');
  if fg[fg_thiophosphoric_acid_halide]     then writeln('thiophosphoric acid halide');
  if fg[fg_thiophosphoric_acid_amide]      then writeln('thiophosphoric acid amide');
  if fg[fg_phosphonic_acid_deriv]          then writeln('phosphonic acid derivative ');
  if fg[fg_phosphonic_acid]                then writeln('phosphonic acid');
  if fg[fg_phosphonic_acid_ester]          then writeln('phosphonic acid ester');
  if fg[fg_phosphine]                      then writeln('phosphine');
  if fg[fg_phosphinoxide]                  then writeln('phosphine oxide');
  if fg[fg_boronic_acid_deriv]             then writeln('boronic acid derivative');
  if fg[fg_boronic_acid]                   then writeln('boronic acid');
  if fg[fg_boronic_acid_ester]             then writeln('boronic acid ester');
  if fg[fg_alkene]                         then writeln('alkene');
  if fg[fg_alkyne]                         then writeln('alkyne');
  if fg[fg_aromatic]                       then writeln('aromatic compound');
  if fg[fg_heterocycle]                    then writeln('heterocyclic compound');
  if fg[fg_alpha_aminoacid]                then writeln('alpha-aminoacid');
  if fg[fg_alpha_hydroxyacid]              then writeln('alpha-hydroxyacid');
end;
*)


procedure write_fg_text(lang:integer);
var
  i : integer;
begin
  for i := 1 to used_fg do
    begin   // first define some exceptions
      if (fg[i] = true) then
        begin
          if (
               (i = fg_carbonyl) or 
               (i = fg_thiocarbonyl) or 
               (i = fg_alcohol) or
               (i = fg_hydroxy) or 
               (i = fg_ether) or 
               (i = fg_amine) or 
               (i = fg_halogen_deriv) or
               (i = fg_alkyl_halide) or 
               (i = fg_aryl_halide) or 
               (i = fg_carboxylic_acid_deriv) or
               (i = fg_carboxylic_acid_amide) or 
               (i = fg_acyl_halide) or 
               (i = fg_thiocarboxylic_acid_deriv) or
               (i = fg_carboxylic_acid_imide) or 
               (i = fg_carbonic_acid_deriv) or 
               (i = fg_carbamic_acid_deriv) or
               (i = fg_thiocarbamic_acid_deriv) or 
               (i = fg_sulfuric_acid_deriv) or 
               (i = fg_sulfonic_acid_deriv) or
               (i = fg_sulfinic_acid_deriv) or 
               (i = fg_sulfenic_acid_deriv) or 
               (i = fg_phosphoric_acid_deriv) or
               (i = fg_thiophosphoric_acid_deriv)
             ) then
            begin
              if (i = fg_hydroxy) and hydroxy_generic    then writeln(mkfglabel(fg_hydroxy,lang));
              if (i = fg_ether) and ether_generic        then writeln(mkfglabel(fg_ether,lang));
              if (i = fg_amine) and amine_generic        then writeln(mkfglabel(fg_amine,lang));
              if (i = fg_halogen_deriv)                 then 
                begin
                  if (not fg[fg_alkyl_halide]) and (not fg[fg_aryl_halide]) and (not fg[fg_acyl_halide]) then
                  writeln(mkfglabel(fg_halogen_deriv,lang));
                end;
              if (i = fg_carbonic_acid_deriv) and not
                (fg[fg_carbonic_acid_monoester] or 
                 fg[fg_carbonic_acid_diester] or
                 fg[fg_carbonic_acid_ester_halide])    then writeln(mkfglabel(fg_carbonic_acid_deriv,lang));
              if (i = fg_carbamic_acid_deriv) and not
                (fg[fg_carbamic_acid] or
                 fg[fg_carbamic_acid_ester] or
                 fg[fg_carbamic_acid_halide])          then writeln(mkfglabel(fg_carbamic_acid_deriv,lang));
              if (i = fg_thiocarbamic_acid_deriv) and not
                 (fg[fg_thiocarbamic_acid] or
                  fg[fg_thiocarbamic_acid_ester] or
                  fg[fg_thiocarbamic_acid_halide])     then writeln(mkfglabel(fg_thiocarbamic_acid_deriv,lang));
            end else writeln(mkfglabel(i,lang));  // now treat the rest normally
        end;  // if fg[i]...
    end;  // for i
end;

(*
procedure write_fg_text_de;   // old version, no longer needed
begin
  if fg[fg_cation]                         then writeln('Kation');
  if fg[fg_anion]                          then writeln('Anion');
//  if fg[fg_carbonyl]                       then writeln('Carbonylverbindung');
  if fg[fg_aldehyde]                       then writeln('Aldehyd');
  if fg[fg_ketone]                         then writeln('Keton');
//  if fg[fg_thiocarbonyl]                   then writeln('Thiocarbonylverbindung');
  if fg[fg_thioaldehyde]                   then writeln('Thioaldehyd');
  if fg[fg_thioketone]                     then writeln('Thioketon');
  if fg[fg_imine]                          then writeln('Imin');
  if fg[fg_hydrazone]                      then writeln('Hydrazon');
  if fg[fg_semicarbazone]                  then writeln('Semicarbazon');
  if fg[fg_thiosemicarbazone]              then writeln('Thiosemicarbazon');
  if fg[fg_oxime]                          then writeln('Oxim');
  if fg[fg_oxime_ether]                    then writeln('Oximether');
  if fg[fg_ketene]                         then writeln('Keten');
  if fg[fg_ketene_acetal_deriv]            then writeln('Keten-Acetal oder Derivat');
  if fg[fg_carbonyl_hydrate]               then writeln('Carbonyl-Hydrat');
  if fg[fg_hemiacetal]                     then writeln('Halbacetal');
  if fg[fg_acetal]                         then writeln('Acetal');
  if fg[fg_hemiaminal]                     then writeln('Halbaminal');
  if fg[fg_aminal]                         then writeln('Aminal');
  if fg[fg_thiohemiaminal]                 then writeln('Thiohalbaminal');
  if fg[fg_thioacetal]                     then writeln('Thioacetal');
  if fg[fg_enamine]                        then writeln('Enamin');
  if fg[fg_enol]                           then writeln('Enol');
  if fg[fg_enolether]                      then writeln('Enolether');
  if fg[fg_hydroxy] and hydroxy_generic    then writeln('Hydroxy-Verbindung');
//  if fg[fg_alcohol]                        then writeln('Alkohol');
  if fg[fg_prim_alcohol]                   then writeln('primärer Alkohol');
  if fg[fg_sec_alcohol]                    then writeln('sekundärer Alkohol');
  if fg[fg_tert_alcohol]                   then writeln('tertiärer Alkohol');
  if fg[fg_1_2_diol]                       then writeln('1,2-Diol');
  if fg[fg_1_2_aminoalcohol]               then writeln('1,2-Aminoalkohol');
  if fg[fg_phenol]                         then writeln('Phenol oder Hydroxyhetaren');
  if fg[fg_1_2_diphenol]                   then writeln('1,2-Diphenol');
  if fg[fg_enediol]                        then writeln('Endiol');
  if fg[fg_ether] and ether_generic        then writeln('Ether');
  if fg[fg_dialkylether]                   then writeln('Dialkylether');
  if fg[fg_alkylarylether]                 then writeln('Alkylarylether ');
  if fg[fg_diarylether]                    then writeln('Diarylether');
  if fg[fg_thioether]                      then writeln('Thioether');
  if fg[fg_disulfide]                      then writeln('Disulfid');
  if fg[fg_peroxide]                       then writeln('Peroxid');
  if fg[fg_hydroperoxide]                  then writeln('Hydroperoxid');
  if fg[fg_hydrazine]                      then writeln('Hydrazin-Derivat');
  if fg[fg_hydroxylamine]                  then writeln('Hydroxylamin');
  if fg[fg_amine] and amine_generic        then writeln('Amin');
  if fg[fg_prim_amine]                     then writeln('primäres Amin');
  if fg[fg_prim_aliph_amine]               then writeln('primäres aliphatisches Amin (Alkylamin)');
  if fg[fg_prim_arom_amine]                then writeln('primäres aromatisches Amin');
  if fg[fg_sec_amine]                      then writeln('sekundäres Amin');
  if fg[fg_sec_aliph_amine]                then writeln('sekundäres aliphatisches Amin (Dialkylamin)');
  if fg[fg_sec_mixed_amine]                then writeln('sekundäres aliphatisches/aromatisches Amin (Alkylarylamin)');
  if fg[fg_sec_arom_amine]                 then writeln('sekundäres aromatisches Amin (Diarylamin)');
  if fg[fg_tert_amine]                     then writeln('tertiäres Amin');
  if fg[fg_tert_aliph_amine]               then writeln('tertiäres aliphatisches Amin (Trialkylamin)');
  if fg[fg_tert_mixed_amine]               then writeln('tertiäres aliphatisches/aromatisches Amin (Alkylarylamin)');
  if fg[fg_tert_arom_amine]                then writeln('tertiäres aromatisches Amin (Triarylamin)');
  if fg[fg_quart_ammonium]                 then writeln('quartäres Ammoniumsalz');
  if fg[fg_n_oxide]                        then writeln('N-Oxid');
  // new in v0.2f
  if fg[fg_halogen_deriv]                  then 
    begin
      if (not fg[fg_alkyl_halide]) and (not fg[fg_aryl_halide]) and (not fg[fg_acyl_halide]) then
      writeln('Halogenverbindung');
    end;
//  if fg[fg_alkyl_halide]                   then writeln('Alkylhalogenid');
  if fg[fg_alkyl_fluoride]                 then writeln('Alkylfluorid');
  if fg[fg_alkyl_chloride]                 then writeln('Alkylchlorid');
  if fg[fg_alkyl_bromide]                  then writeln('Alkylbromid');
  if fg[fg_alkyl_iodide]                   then writeln('Alkyliodid');
//  if fg[fg_aryl_halide]                    then writeln('Arylhalogenid');
  if fg[fg_aryl_fluoride]                  then writeln('Arylfluorid');
  if fg[fg_aryl_chloride]                  then writeln('Arylchlorid');
  if fg[fg_aryl_bromide]                   then writeln('Arylbromid');
  if fg[fg_aryl_iodide]                    then writeln('Aryliodid');
  if fg[fg_organometallic]                 then writeln('Organometall-Verbindung');
  if fg[fg_organolithium]                  then writeln('Organolithium-Verbindung');
  if fg[fg_organomagnesium]                then writeln('Organomagnesium-Verbindung');
//  if fg[fg_carboxylic_acid_deriv]          then writeln('Carbonsäure-Derivat');
  if fg[fg_carboxylic_acid]                then writeln('Carbonsäure');
  if fg[fg_carboxylic_acid_salt]           then writeln('Carbonsäuresalz');
  if fg[fg_carboxylic_acid_ester]          then writeln('Carbonsäureester');
  if fg[fg_lactone]                        then writeln('Lacton');
//  if fg[fg_carboxylic_acid_amide]          then writeln('Carbonsäureamid');
  if fg[fg_carboxylic_acid_prim_amide]     then writeln('primäres Carbonsäureamid');
  if fg[fg_carboxylic_acid_sec_amide]      then writeln('sekundäres Carbonsäureamid');
  if fg[fg_carboxylic_acid_tert_amide]     then writeln('tertiäres Carbonsäureamid');
  if fg[fg_lactam]                         then writeln('Lactam');
  if fg[fg_carboxylic_acid_hydrazide]      then writeln('Carbonsäurehydrazid');
  if fg[fg_carboxylic_acid_azide]          then writeln('Carbonsäureazid');
  if fg[fg_hydroxamic_acid]                then writeln('Hydroxamsäure');
  if fg[fg_carboxylic_acid_amidine]        then writeln('Carbonsäureamidin');
  if fg[fg_carboxylic_acid_amidrazone]     then writeln('Carbonsäureamidrazon');
  if fg[fg_nitrile]                        then writeln('Carbonitril');
//  if fg[fg_acyl_halide]                    then writeln('Acylhalogenid');
  if fg[fg_acyl_fluoride]                  then writeln('Acylfluorid');
  if fg[fg_acyl_chloride]                  then writeln('Acylchlorid');
  if fg[fg_acyl_bromide]                   then writeln('Acylbromid');
  if fg[fg_acyl_iodide]                    then writeln('Acyliodid');
  if fg[fg_acyl_cyanide]                   then writeln('Acylcyanid');
  if fg[fg_imido_ester]                    then writeln('Imidoester');
  if fg[fg_imidoyl_halide]                 then writeln('Imidoylhalogenid');
//  if fg[fg_thiocarboxylic_acid_deriv]      then writeln('Thiocarbonsäure-Derivat');
  if fg[fg_thiocarboxylic_acid]            then writeln('Thiocarbonsäure');
  if fg[fg_thiocarboxylic_acid_ester]      then writeln('Thiocarbonsäureester');
  if fg[fg_thiolactone]                    then writeln('Thiolacton');
  if fg[fg_thiocarboxylic_acid_amide]      then writeln('Thiocarbonsäureamid');
  if fg[fg_thiolactam]                     then writeln('Thiolactam');
  if fg[fg_imido_thioester]                then writeln('Imidothioester');
  if fg[fg_oxohetarene]                    then writeln('Oxo(het)aren');
  if fg[fg_thioxohetarene]                 then writeln('Thioxo(het)aren');
  if fg[fg_iminohetarene]                  then writeln('Imino(het)aren');
  if fg[fg_orthocarboxylic_acid_deriv]     then writeln('Orthocarbonsäure-Derivat');
  if fg[fg_carboxylic_acid_orthoester]     then writeln('Orthoester');
  if fg[fg_carboxylic_acid_amide_acetal]   then writeln('Amidacetal');
  if fg[fg_carboxylic_acid_anhydride]      then writeln('Carbonsäureanhydrid');
//  if fg[fg_carboxylic_acid_imide]          then writeln('Carbonsäureimid');
  if fg[fg_carboxylic_acid_unsubst_imide]  then writeln('Carbonsäureimid, N-unsubstituiert');
  if fg[fg_carboxylic_acid_subst_imide]    then writeln('Carbonsäureimid, N-substituiert');
  if fg[fg_co2_deriv]                      then writeln('CO2-Derivat (allgemein)');
  if fg[fg_carbonic_acid_deriv] and not    // changed in v0.3c
    (fg[fg_carbonic_acid_monoester] or 
     fg[fg_carbonic_acid_diester] or
     fg[fg_carbonic_acid_ester_halide])    then writeln('Kohlensäure-Derivat');
  if fg[fg_carbonic_acid_monoester]        then writeln('Kohlensäuremonoester');
  if fg[fg_carbonic_acid_diester]          then writeln('Kohlensäurediester');
  if fg[fg_carbonic_acid_ester_halide]     then writeln('Kohlensäureesterhalogenid (Alkyl/Aryl-Halogenformiat)');
  if fg[fg_thiocarbonic_acid_deriv]        then writeln('Thiokohlensäure-Derivat');
  if fg[fg_thiocarbonic_acid_monoester]    then writeln('Thiokohlensäuremonoester');
  if fg[fg_thiocarbonic_acid_diester]      then writeln('Thiokohlensäurediester');
  if fg[fg_thiocarbonic_acid_ester_halide] then writeln('Thiokohlensäureesterhalogenid (Alkyl/Aryl-Halogenthioformiat)');
  if fg[fg_carbamic_acid_deriv] and not    // changed in v0.3c
    (fg[fg_carbamic_acid] or
     fg[fg_carbamic_acid_ester] or
     fg[fg_carbamic_acid_halide])          then writeln('Carbaminsäure-Derivat');
  if fg[fg_carbamic_acid]                  then writeln('Carbaminsäure');
  if fg[fg_carbamic_acid_ester]            then writeln('Carbaminsäureester (Urethan)');
  if fg[fg_carbamic_acid_halide]           then writeln('Carbaminsäurehalogenid (Halogenformamid)');
  if fg[fg_thiocarbamic_acid_deriv] and not  // changed in v0.3c
     (fg[fg_thiocarbamic_acid] or
      fg[fg_thiocarbamic_acid_ester] or
      fg[fg_thiocarbamic_acid_halide])     then writeln('Thiocarbaminsäure-Derivat');
  if fg[fg_thiocarbamic_acid]              then writeln('Thiocarbaminsäure');
  if fg[fg_thiocarbamic_acid_ester]        then writeln('Thiocarbaminsäureester');
  if fg[fg_thiocarbamic_acid_halide]       then writeln('Thiocarbaminsäurehalogenid (Halogenthioformamid)');
  if fg[fg_urea]                           then writeln('Harnstoff');
  if fg[fg_isourea]                        then writeln('Isoharnstoff');
  if fg[fg_thiourea]                       then writeln('Thioharnstoff');
  if fg[fg_isothiourea]                    then writeln('Isothioharnstoff');
  if fg[fg_guanidine]                      then writeln('Guanidin');
  if fg[fg_semicarbazide]                  then writeln('Semicarbazid');
  if fg[fg_thiosemicarbazide]              then writeln('Thiosemicarbazid');
  if fg[fg_azide]                          then writeln('Azid');
  if fg[fg_azo_compound]                   then writeln('Azoverbindung');
  if fg[fg_diazonium_salt]                 then writeln('Diazoniumsalz');
  if fg[fg_isonitrile]                     then writeln('Isonitril');
  if fg[fg_cyanate]                        then writeln('Cyanat');
  if fg[fg_isocyanate]                     then writeln('Isocyanat');
  if fg[fg_thiocyanate]                    then writeln('Thiocyanat');
  if fg[fg_isothiocyanate]                 then writeln('Isothiocyanat');
  if fg[fg_carbodiimide]                   then writeln('Carbodiimid');
  if fg[fg_nitroso_compound]               then writeln('Nitroso-Verbindung');
  if fg[fg_nitro_compound]                 then writeln('Nitro-Verbindung');
  if fg[fg_nitrite]                        then writeln('Nitrit');
  if fg[fg_nitrate]                        then writeln('Nitrat');
//  if fg[fg_sulfuric_acid_deriv]            then writeln('Schwefelsäure-Derivat');
  if fg[fg_sulfuric_acid]                  then writeln('Schwefelsäure');
  if fg[fg_sulfuric_acid_monoester]        then writeln('Schwefelsäuremonoester');
  if fg[fg_sulfuric_acid_diester]          then writeln('Schwefelsäurediester');
  if fg[fg_sulfuric_acid_amide_ester]      then writeln('Schwefelsäureamidester');
  if fg[fg_sulfuric_acid_amide]            then writeln('Schwefelsäureamid');
  if fg[fg_sulfuric_acid_diamide]          then writeln('Schwefelsäurediamid');
  if fg[fg_sulfuryl_halide]                then writeln('Sulfurylhalogenid');
//  if fg[fg_sulfonic_acid_deriv]            then writeln('Sulfonsäure-Derivat ');
  if fg[fg_sulfonic_acid]                  then writeln('Sulfonsäure');
  if fg[fg_sulfonic_acid_ester]            then writeln('Sulfonsäureester');
  if fg[fg_sulfonamide]                    then writeln('Sulfonamid');
  if fg[fg_sulfonyl_halide]                then writeln('Sulfonylhalogenid');
  if fg[fg_sulfone]                        then writeln('Sulfon');
  if fg[fg_sulfoxide]                      then writeln('Sulfoxid');
//  if fg[fg_sulfinic_acid_deriv]            then writeln('Sulfinsäure-Derivat');
  if fg[fg_sulfinic_acid]                  then writeln('Sulfinsäure');
  if fg[fg_sulfinic_acid_ester]            then writeln('Sulfinsäureester');
  if fg[fg_sulfinic_acid_halide]           then writeln('Sulfinsäurehalogenid');
  if fg[fg_sulfinic_acid_amide]            then writeln('Sulfinsäureamid');
//  if fg[fg_sulfenic_acid_deriv]            then writeln('Sulfensäure-Derivat');
  if fg[fg_sulfenic_acid]                  then writeln('Sulfensäure');
  if fg[fg_sulfenic_acid_ester]            then writeln('Sulfensäureester');
  if fg[fg_sulfenic_acid_halide]           then writeln('Sulfensäurehalogenid');
  if fg[fg_sulfenic_acid_amide]            then writeln('Sulfensäureamid');
  if fg[fg_thiol]                          then writeln('Thiol (Sulfanyl-Verbindung, Mercaptan)');
  if fg[fg_alkylthiol]                     then writeln('Alkylthiol');
  if fg[fg_arylthiol]                      then writeln('Arylthiol');
//  if fg[fg_phosphoric_acid_deriv]          then writeln('Phosphorsäure-Derivat');
  if fg[fg_phosphoric_acid]                then writeln('Phosphorsäure');
  if fg[fg_phosphoric_acid_ester]          then writeln('Phosphorsäureester');
  if fg[fg_phosphoric_acid_halide]         then writeln('Phosphorsäurehalogenid');
  if fg[fg_phosphoric_acid_amide]          then writeln('Phosphorsäureamid');
//  if fg[fg_thiophosphoric_acid_deriv]      then writeln('Thiophosphorsäure-Derivat');
  if fg[fg_thiophosphoric_acid]            then writeln('Thiophosphorsäure');
  if fg[fg_thiophosphoric_acid_ester]      then writeln('Thiophosphorsäureester');
  if fg[fg_thiophosphoric_acid_halide]     then writeln('Thiophosphorsäurehalogenid');
  if fg[fg_thiophosphoric_acid_amide]      then writeln('Thiophosphorsäureamid');
  if fg[fg_phosphonic_acid_deriv]          then writeln('Phosphonsäure-Derivat ');
  if fg[fg_phosphonic_acid]                then writeln('Phosphonsäure');
  if fg[fg_phosphonic_acid_ester]          then writeln('Phosphonsäureester');
  if fg[fg_phosphine]                      then writeln('Phosphin');
  if fg[fg_phosphinoxide]                  then writeln('Phosphinoxid');
  if fg[fg_boronic_acid_deriv]             then writeln('Boronsäure-Derivat');
  if fg[fg_boronic_acid]                   then writeln('Boronsäure');
  if fg[fg_boronic_acid_ester]             then writeln('Boronsäureester');
  if fg[fg_alkene]                         then writeln('Alken');
  if fg[fg_alkyne]                         then writeln('Alkin');
  if fg[fg_aromatic]                       then writeln('aromatische Verbindung');
  if fg[fg_heterocycle]                    then writeln('heterocyclische Verbindung');
  if fg[fg_alpha_aminoacid]                then writeln('alpha-Aminosäure');
  if fg[fg_alpha_hydroxyacid]              then writeln('alpha-Hydroxysäure');
end;
*)

procedure write_fg_code;
const
  sc = ';';
var
  i : integer;
begin
  for i := 1 to used_fg do
    begin   // first define some exceptions
      if (fg[i] = true) then
        begin
          if (
               (i = fg_carbonyl) or 
               (i = fg_thiocarbonyl) or 
               (i = fg_alcohol) or
               (i = fg_hydroxy) or 
               (i = fg_ether) or 
               (i = fg_amine) or 
               (i = fg_prim_amine) or 
               (i = fg_sec_amine) or 
               (i = fg_tert_amine) or 
               (i = fg_halogen_deriv) or
               (i = fg_alkyl_halide) or 
               (i = fg_aryl_halide) or 
               (i = fg_carboxylic_acid_deriv) or
               (i = fg_carboxylic_acid_amide) or 
               (i = fg_acyl_halide) or 
               (i = fg_thiocarboxylic_acid_deriv) or
               (i = fg_carboxylic_acid_imide) 
             ) then
            begin
              if (i = fg_hydroxy) and hydroxy_generic    then write(mkfglabel(fg_hydroxy,0),sc);
              if (i = fg_ether) and ether_generic        then write(mkfglabel(fg_ether,0),sc);
              if (i = fg_amine) and amine_generic        then write(mkfglabel(fg_amine,0),sc);
              if (i = fg_halogen_deriv)                  then 
                begin
                  if (not fg[fg_alkyl_halide]) and (not fg[fg_aryl_halide]) and (not fg[fg_acyl_halide]) then
                  write(mkfglabel(fg_halogen_deriv,0),sc);
                end;
            end else write(mkfglabel(i,0),sc);  // now treat the rest normally
        end;  // if fg[i]...
    end;  // for i
  writeln;
end;


(*
procedure write_fg_code;
const
  sc = ';';
begin
  if fg[fg_cation]                         then write('000000T2',sc);
  if fg[fg_anion]                          then write('000000T1',sc);
//  if fg[fg_carbonyl]                       then write('C2O10000',sc);
  if fg[fg_aldehyde]                       then write('C2O1H000',sc);
  if fg[fg_ketone]                         then write('C2O1C000',sc);
//  if fg[fg_thiocarbonyl]                   then write('C2S10000',sc);
  if fg[fg_thioaldehyde]                   then write('C2S1H000',sc);
  if fg[fg_thioketone]                     then write('C2S1C000',sc);
  if fg[fg_imine]                          then write('C2N10000',sc);
  if fg[fg_hydrazone]                      then write('C2N1N000',sc);
  if fg[fg_semicarbazone]                  then write('C2NNC4ON',sc);
  if fg[fg_thiosemicarbazone]              then write('C2NNC4SN',sc);
  if fg[fg_oxime]                          then write('C2N1OH00',sc);
  if fg[fg_oxime_ether]                    then write('C2N1OC00',sc);
  if fg[fg_ketene]                         then write('C3OC0000',sc);
  if fg[fg_ketene_acetal_deriv]            then write('C3OCC000',sc);
  if fg[fg_carbonyl_hydrate]               then write('C2O2H200',sc);
  if fg[fg_hemiacetal]                     then write('C2O2HC00',sc);
  if fg[fg_acetal]                         then write('C2O2CC00',sc);
  if fg[fg_hemiaminal]                     then write('C2NOHC10',sc);
  if fg[fg_aminal]                         then write('C2N2CC10',sc);
  if fg[fg_thiohemiaminal]                 then write('C2NSHC10',sc);
  if fg[fg_thioacetal]                     then write('C2S2CC00',sc);
  if fg[fg_enamine]                        then write('C2CNH000',sc);
  if fg[fg_enol]                           then write('C2COH000',sc);
  if fg[fg_enolether]                      then write('C2COC000',sc);
  if fg[fg_hydroxy] and hydroxy_generic    then write('O1H00000',sc);
//  if fg[fg_alcohol]                        then write('O1H0C000',sc);
  if fg[fg_prim_alcohol]                   then write('O1H1C000',sc);
  if fg[fg_sec_alcohol]                    then write('O1H2C000',sc);
  if fg[fg_tert_alcohol]                   then write('O1H3C000',sc);
  if fg[fg_1_2_diol]                       then write('O1H0CO1H',sc);
  if fg[fg_1_2_aminoalcohol]               then write('O1H0CN1C',sc);
  if fg[fg_phenol]                         then write('O1H1A000',sc);
  if fg[fg_1_2_diphenol]                   then write('O1H2A000',sc);
  if fg[fg_enediol]                        then write('C2COH200',sc);
  if fg[fg_ether] and ether_generic        then write('O1C00000',sc);
  if fg[fg_dialkylether]                   then write('O1C0CC00',sc);
  if fg[fg_alkylarylether]                 then write('O1C0CA00',sc);
  if fg[fg_diarylether]                    then write('O1C0AA00',sc);
  if fg[fg_thioether]                      then write('S1C00000',sc);
  if fg[fg_disulfide]                      then write('S1S1C000',sc);
  if fg[fg_peroxide]                       then write('O1O1C000',sc);
  if fg[fg_hydroperoxide]                  then write('O1O1H000',sc);
  if fg[fg_hydrazine]                      then write('N1N10000',sc);
  if fg[fg_hydroxylamine]                  then write('N1O1H000',sc);
  if fg[fg_amine] and amine_generic        then write('N1C00000',sc);
//  if fg[fg_prim_amine]                     then write('N1C10000',sc);
  if fg[fg_prim_aliph_amine]               then write('N1C1C000',sc);
  if fg[fg_prim_arom_amine]                then write('N1C1A000',sc);
//  if fg[fg_sec_amine]                      then write('N1C20000',sc);
  if fg[fg_sec_aliph_amine]                then write('N1C2CC00',sc);
  if fg[fg_sec_mixed_amine]                then write('N1C2AC00',sc);
  if fg[fg_sec_arom_amine]                 then write('N1C2AA00',sc);
//  if fg[fg_tert_amine]                     then write('N1C30000',sc);
  if fg[fg_tert_aliph_amine]               then write('N1C3CC00',sc);
  if fg[fg_tert_mixed_amine]               then write('N1C3AC00',sc);
  if fg[fg_tert_arom_amine]                then write('N1C3AA00',sc);
  if fg[fg_quart_ammonium]                 then write('N1C400T2',sc);
  if fg[fg_n_oxide]                        then write('N0O10000',sc);
//  if fg[fg_halogen_deriv]                  then write('XX000000',sc);
  // new in v0.2f
  if fg[fg_halogen_deriv]                  then 
    begin
      if (not fg[fg_alkyl_halide]) and (not fg[fg_aryl_halide]) and (not fg[fg_acyl_halide]) then
      write('XX000000',sc);
    end;
//  if fg[fg_alkyl_halide]                   then write('XX00C000',sc);
  if fg[fg_alkyl_fluoride]                 then write('XF00C000',sc);
  if fg[fg_alkyl_chloride]                 then write('XC00C000',sc);
  if fg[fg_alkyl_bromide]                  then write('XB00C000',sc);
  if fg[fg_alkyl_iodide]                   then write('XI00C000',sc);
//  if fg[fg_aryl_halide]                    then write('XX00A000',sc);
  if fg[fg_aryl_fluoride]                  then write('XF00A000',sc);
  if fg[fg_aryl_chloride]                  then write('XC00A000',sc);
  if fg[fg_aryl_bromide]                   then write('XB00A000',sc);
  if fg[fg_aryl_iodide]                    then write('XI00A000',sc);
  if fg[fg_organometallic]                 then write('000000MX',sc);
  if fg[fg_organolithium]                  then write('000000ML',sc);
  if fg[fg_organomagnesium]                then write('000000MM',sc);
//  if fg[fg_carboxylic_acid_deriv]          then write('C3O20000',sc);
  if fg[fg_carboxylic_acid]                then write('C3O2H000',sc);
  if fg[fg_carboxylic_acid_salt]           then write('C3O200T1',sc);
  if fg[fg_carboxylic_acid_ester]          then write('C3O2C000',sc);
  if fg[fg_lactone]                        then write('C3O2CZ00',sc);
//  if fg[fg_carboxylic_acid_amide]          then write('C3ONC000',sc);
  if fg[fg_carboxylic_acid_prim_amide]     then write('C3ONC100',sc);
  if fg[fg_carboxylic_acid_sec_amide]      then write('C3ONC200',sc);
  if fg[fg_carboxylic_acid_tert_amide]     then write('C3ONC300',sc);
  if fg[fg_lactam]                         then write('C3ONCZ00',sc);
  if fg[fg_carboxylic_acid_hydrazide]      then write('C3ONN100',sc);
  if fg[fg_carboxylic_acid_azide]          then write('C3ONN200',sc);
  if fg[fg_hydroxamic_acid]                then write('C3ONOH00',sc);
  if fg[fg_carboxylic_acid_amidine]        then write('C3N2H000',sc);
  if fg[fg_carboxylic_acid_amidrazone]     then write('C3NNN100',sc);
  if fg[fg_nitrile]                        then write('C3N00000',sc);
//  if fg[fg_acyl_halide]                    then write('C3OXX000',sc);
  if fg[fg_acyl_fluoride]                  then write('C3OXF000',sc);
  if fg[fg_acyl_chloride]                  then write('C3OXC000',sc);
  if fg[fg_acyl_bromide]                   then write('C3OXB000',sc);
  if fg[fg_acyl_iodide]                    then write('C3OXI000',sc);
  if fg[fg_acyl_cyanide]                   then write('C2OC3N00',sc);
  if fg[fg_imido_ester]                    then write('C3NOC000',sc);
  if fg[fg_imidoyl_halide]                 then write('C3NXX000',sc);
//  if fg[fg_thiocarboxylic_acid_deriv]      then write('C3SO0000',sc);
  if fg[fg_thiocarboxylic_acid]            then write('C3SOH000',sc);
  if fg[fg_thiocarboxylic_acid_ester]      then write('C3SOC000',sc);
  if fg[fg_thiolactone]                    then write('C3SOCZ00',sc);
  if fg[fg_thiocarboxylic_acid_amide]      then write('C3SNH000',sc);
  if fg[fg_thiolactam]                     then write('C3SNCZ00',sc);
  if fg[fg_imido_thioester]                then write('C3NSC000',sc);
  if fg[fg_oxohetarene]                    then write('C3ONAZ00',sc);
  if fg[fg_thioxohetarene]                 then write('C3SNAZ00',sc);
  if fg[fg_iminohetarene]                  then write('C3NNAZ00',sc);
  if fg[fg_orthocarboxylic_acid_deriv]     then write('C3O30000',sc);
  if fg[fg_carboxylic_acid_orthoester]     then write('C3O3C000',sc);
  if fg[fg_carboxylic_acid_amide_acetal]   then write('C3O3NC00',sc);
  if fg[fg_carboxylic_acid_anhydride]      then write('C3O2C3O2',sc);
//  if fg[fg_carboxylic_acid_imide]          then write('C3ONC000',sc);
  if fg[fg_carboxylic_acid_unsubst_imide]  then write('C3ONCH10',sc);
  if fg[fg_carboxylic_acid_subst_imide]    then write('C3ONCC10',sc);
  if fg[fg_co2_deriv]                      then write('C4000000',sc);
  if fg[fg_carbonic_acid_deriv]            then write('C4O30000',sc);
  if fg[fg_carbonic_acid_monoester]        then write('C4O3C100',sc);
  if fg[fg_carbonic_acid_diester]          then write('C4O3C200',sc);
  if fg[fg_carbonic_acid_ester_halide]     then write('C4O3CX00',sc);
  if fg[fg_thiocarbonic_acid_deriv]        then write('C4SO0000',sc);
  if fg[fg_thiocarbonic_acid_monoester]    then write('C4SOC100',sc);
  if fg[fg_thiocarbonic_acid_diester]      then write('C4SOC200',sc);
  if fg[fg_thiocarbonic_acid_ester_halide] then write('C4SOX_00',sc);
  if fg[fg_carbamic_acid_deriv]            then write('C4O2N000',sc);
  if fg[fg_carbamic_acid]                  then write('C4O2NH00',sc);
  if fg[fg_carbamic_acid_ester]            then write('C4O2NC00',sc);
  if fg[fg_carbamic_acid_halide]           then write('C4O2NX00',sc);
  if fg[fg_thiocarbamic_acid_deriv]        then write('C4SN0000',sc);
  if fg[fg_thiocarbamic_acid]              then write('C4SNOH00',sc);
  if fg[fg_thiocarbamic_acid_ester]        then write('C4SNOC00',sc);
  if fg[fg_thiocarbamic_acid_halide]       then write('C4SNXX00',sc);
  if fg[fg_urea]                           then write('C4O1N200',sc);
  if fg[fg_isourea]                        then write('C4N2O100',sc);
  if fg[fg_thiourea]                       then write('C4S1N200',sc);
  if fg[fg_isothiourea]                    then write('C4N2S100',sc);
  if fg[fg_guanidine]                      then write('C4N30000',sc);
  if fg[fg_semicarbazide]                  then write('C4ON2N00',sc);
  if fg[fg_thiosemicarbazide]              then write('C4SN2N00',sc);
  if fg[fg_azide]                          then write('N4N20000',sc);
  if fg[fg_azo_compound]                   then write('N2N10000',sc);
  if fg[fg_diazonium_salt]                 then write('N3N100T2',sc);
  if fg[fg_isonitrile]                     then write('N3C10000',sc);
  if fg[fg_cyanate]                        then write('C4NO1000',sc);
  if fg[fg_isocyanate]                     then write('C4NO2000',sc);
  if fg[fg_thiocyanate]                    then write('C4NS1000',sc);
  if fg[fg_isothiocyanate]                 then write('C4NS2000',sc);
  if fg[fg_carbodiimide]                   then write('C4N20000',sc);
  if fg[fg_nitroso_compound]               then write('N2O10000',sc);
  if fg[fg_nitro_compound]                 then write('N4O20000',sc);
  if fg[fg_nitrite]                        then write('N3O20000',sc);
  if fg[fg_nitrate]                        then write('N4O30000',sc);
  if fg[fg_sulfuric_acid_deriv]            then write('S6O00000',sc);
  if fg[fg_sulfuric_acid]                  then write('S6O4H000',sc);
  if fg[fg_sulfuric_acid_monoester]        then write('S6O4HC00',sc);
  if fg[fg_sulfuric_acid_diester]          then write('S6O4CC00',sc);
  if fg[fg_sulfuric_acid_amide_ester]      then write('S6O3NC00',sc);
  if fg[fg_sulfuric_acid_amide]            then write('S6O3N100',sc);
  if fg[fg_sulfuric_acid_diamide]          then write('S6O2N200',sc);
  if fg[fg_sulfuryl_halide]                then write('S6O3XX00',sc);
  if fg[fg_sulfonic_acid_deriv]            then write('S5O00000',sc);
  if fg[fg_sulfonic_acid]                  then write('S5O3H000',sc);
  if fg[fg_sulfonic_acid_ester]            then write('S5O3C000',sc);
  if fg[fg_sulfonamide]                    then write('S5O2N000',sc);
  if fg[fg_sulfonyl_halide]                then write('S5O2XX00',sc);
  if fg[fg_sulfone]                        then write('S4O20000',sc);
  if fg[fg_sulfoxide]                      then write('S2O10000',sc);
  if fg[fg_sulfinic_acid_deriv]            then write('S3O00000',sc);
  if fg[fg_sulfinic_acid]                  then write('S3O2H000',sc);
  if fg[fg_sulfinic_acid_ester]            then write('S3O2C000',sc);
  if fg[fg_sulfinic_acid_halide]           then write('S3O1XX00',sc);
  if fg[fg_sulfinic_acid_amide]            then write('S3O1N000',sc);
  if fg[fg_sulfenic_acid_deriv]            then write('S1O00000',sc);
  if fg[fg_sulfenic_acid]                  then write('S1O1H000',sc);
  if fg[fg_sulfenic_acid_ester]            then write('S1O1C000',sc);
  if fg[fg_sulfenic_acid_halide]           then write('S1O0XX00',sc);
  if fg[fg_sulfenic_acid_amide]            then write('S1O0N100',sc);
//  if fg[fg_thiol]                          then write('S1H10000',sc);
  if fg[fg_alkylthiol]                     then write('S1H1C000',sc);
  if fg[fg_arylthiol]                      then write('S1H1A000',sc);
  if fg[fg_phosphoric_acid_deriv]          then write('P5O0H000',sc);
  if fg[fg_phosphoric_acid]                then write('P5O4H200',sc);
  if fg[fg_phosphoric_acid_ester]          then write('P5O4HC00',sc);
  if fg[fg_phosphoric_acid_halide]         then write('P5O3HX00',sc);
  if fg[fg_phosphoric_acid_amide]          then write('P5O3HN00',sc);
  if fg[fg_thiophosphoric_acid_deriv]      then write('P5O0S000',sc);
  if fg[fg_thiophosphoric_acid]            then write('P5O3SH00',sc);
  if fg[fg_thiophosphoric_acid_ester]      then write('P5O3SC00',sc);
  if fg[fg_thiophosphoric_acid_halide]     then write('P5O2SX00',sc);
  if fg[fg_thiophosphoric_acid_amide]      then write('P5O2SN00',sc);
  if fg[fg_phosphonic_acid_deriv]          then write('P4O30000',sc);
  if fg[fg_phosphonic_acid]                then write('P4O3H000',sc);
  if fg[fg_phosphonic_acid_ester]          then write('P4O3C000',sc);
  if fg[fg_phosphine]                      then write('P3000000',sc);
  if fg[fg_phosphinoxide]                  then write('P2O00000',sc);
  if fg[fg_boronic_acid_deriv]             then write('B2O20000',sc);
  if fg[fg_boronic_acid]                   then write('B2O2H000',sc);
  if fg[fg_boronic_acid_ester]             then write('B2O2C000',sc);
  if fg[fg_alkene]                         then write('000C2C00',sc);
  if fg[fg_alkyne]                         then write('000C3C00',sc);
  if fg[fg_aromatic]                       then write('0000A000',sc);
  if fg[fg_heterocycle]                    then write('0000CZ00',sc);
  if fg[fg_alpha_aminoacid]                then write('C3O2HN1C',sc);
  if fg[fg_alpha_hydroxyacid]              then write('C3O2HO1H',sc);
  writeln;  // v0.4
end;
*)

(*
procedure write_fg_binary;
var
  i : integer;
  n : integer;
  o : char;
begin
  for i := 1 to (max_fg div 8) do
    begin
      n := 0;
      if fg[8*i]    then inc(n);
      if fg[8*i-1]  then inc(n,2);
      if fg[8*i-2]  then inc(n,4);
      if fg[8*i-3]  then inc(n,8);
      if fg[8*i-4]  then inc(n,16);
      if fg[8*i-5]  then inc(n,32);
      if fg[8*i-6]  then inc(n,64);
      if fg[8*i-7]  then inc(n,128);
      o := chr(n);
      write(o);
    end;
end;
*)

procedure write_fg_binary;  // v0.4
const
  bsize = 32;
var
  i, j, n1 : integer;
  fgincrement : int64;
  fgdecimal : int64;
begin
  n1 := 0;
  for i := 0 to ((max_fg div bsize) - 1) do
    begin
      fgdecimal := 0;
      for j := 1 to bsize do
        begin
          fgincrement := 1;
          if fg[((bsize*i)+j)] then 
            begin
              inc(n1);
              fgincrement := fgincrement shl (j-1);
              fgdecimal := fgdecimal + fgincrement;
            end;
        end;
      if (i > 0) then write(',');
      write(fgdecimal);
    end;
  write(';',n1);
  writeln;
end;


procedure write_fg_bitstring;
var
  i : integer;
begin
  for i := 1 to max_fg do if fg[i] then write('1') else write('0');
  writeln;  // v0.4
end;


procedure write_fg_pos;
var
  i,j,k,n,b_id,r_id,s : integer;
  i_str : string;
begin
  for i := 1 to used_fg do
    begin
      n := fgloc^[i,0];
      if (n > 0) then
        begin
          i_str := inttostr(i);
          while length(i_str) < 3 do i_str := '0' + i_str;
          i_str := '#' + i_str + ':';
          write(i_str);
          if (fglang > 0) then write(mkfglabel(i,fglang),':');
          write(n,':');
          for j := 1 to n do
            begin
              if (i = fg_aromatic ) or (i = fg_heterocycle) then
                begin  // use ring list
                  r_id := fgloc^[i,j];
                  s := ringprop^[r_id].size;
                  if (j > 1) then write(',');
                  for k := 1 to s do
                    begin
                      write(ring^[r_id,k],'-');
                    end;          
                end else
                begin
                  if (i = fg_1_2_diol) or
                     (i = fg_1_2_diphenol) or
                     (i = fg_enediol) or
                     (i = fg_disulfide) or
                     (i = fg_peroxide) or
                     (i = fg_hydrazine) or
                     (i = fg_azo_compound) or
                     (i = fg_alkene) or
                     (i = fg_alkyne) then
                    begin  // use bond list
                      b_id := fgloc^[i,j];
                      if (j > 1) then write(',');
                      write(bond^[b_id].a1,'-',bond^[b_id].a2);
                    end else
                    begin  // use atom list
                      if (j > 1) then write(',');
                      write(fgloc^[i,j]);
                    end;
                end;
            end;
          writeln;
        end;
    end;
end;


procedure readinputfile(molfilename:string);  // new version in v0.2g, updated in v0.4b
var
  rline : string;
begin
  molbufindex := 0;
  if not opt_stdin then
    begin
      if not rfile_is_open then
        begin
          assign(rfile,molfilename);
          reset(rfile);
          rfile_is_open := true;
        end;
      rline := '';
      mol_in_queue := false;
      while (not eof(rfile)) and (pos('$$$$',rline) = 0) do
        begin
          readln(rfile,rline);
          //mol_in_queue := false;
          if molbufindex < (max_atoms+max_bonds+64) then
            begin
              inc(molbufindex);
              molbuf^[molbufindex] := rline;
            end else
            begin
              writeln('Not enough memory for molfile! ',molbufindex);
              close(rfile);
              halt(1);
            end;
          if pos('$$$$',rline) > 0 then 
            begin
              mol_in_queue := true;
              sep_label := get_sep_label(rline);  // v0.4b
            end;
        end;
      if eof(rfile) then
        begin
          close(rfile);
          rfile_is_open := false;
          mol_in_queue := false;
        end;
    end else              // read from standard input
    begin
      rline := '';
      mol_in_queue := false;
      while (not eof) and (pos('$$$$',rline) = 0) do
        begin
          readln(rline);
          if molbufindex < (max_atoms+max_bonds+64) then
            begin
              inc(molbufindex);
              molbuf^[molbufindex] := rline;
            end else
            begin
              writeln('Not enough memory!');
              halt(1);
            end;
          if pos('$$$$',rline) > 0 then 
            begin
              mol_in_queue := true;
              sep_label := get_sep_label(rline);  // v0.4b
            end;
        end;
    end;
end;

procedure copy_mol_to_needle;
var
  i, j : integer;
begin
  if (n_atoms = 0) then exit;
  try
    getmem(ndl_atom,n_atoms*sizeof(atom_rec));
    getmem(ndl_bond,n_bonds*sizeof(bond_rec));
    getmem(ndl_ring,sizeof(ringlist));
    getmem(ndl_ringprop,sizeof(ringprop_type));
  except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end;
  ndl_n_atoms := n_atoms;
  ndl_n_bonds := n_bonds;
  ndl_n_rings := n_rings;
  ndl_n_heavyatoms := n_heavyatoms;
  ndl_n_trueheavyatoms := n_trueheavyatoms;  // v0.4b
  ndl_n_heavybonds := n_heavybonds;
  ndl_molname := molname;
  ndl_n_Ctot := n_Ctot;
  ndl_n_Otot := n_Otot;
  ndl_n_Ntot := n_Ntot;
  for i := 1 to n_atoms do
    begin
      ndl_atom^[i].element        := atom^[i].element;
      ndl_atom^[i].atype          := atom^[i].atype;
      ndl_atom^[i].x              := atom^[i].x;
      ndl_atom^[i].y              := atom^[i].y;
      ndl_atom^[i].z              := atom^[i].z;
      ndl_atom^[i].formal_charge  := atom^[i].formal_charge;
      ndl_atom^[i].real_charge    := atom^[i].real_charge;
      ndl_atom^[i].Hexp           := atom^[i].Hexp;
      ndl_atom^[i].Htot           := atom^[i].Htot;
      ndl_atom^[i].neighbor_count := atom^[i].neighbor_count;
      ndl_atom^[i].ring_count     := atom^[i].ring_count;
      ndl_atom^[i].arom           := atom^[i].arom;
      ndl_atom^[i].stereo_care    := atom^[i].stereo_care;
      ndl_atom^[i].heavy          := atom^[i].heavy;  // v0.3l
      ndl_atom^[i].metal          := atom^[i].metal;  // v0.3l
      ndl_atom^[i].tag            := atom^[i].tag;  // v0.3o
      ndl_atom^[i].nucleon_number := atom^[i].nucleon_number;  // v0.3p
      ndl_atom^[i].radical_type   := atom^[i].radical_type;  // v0.3p
    end;
  if (n_bonds > 0) then
    begin
      for i := 1 to n_bonds do
        begin
          ndl_bond^[i].a1         := bond^[i].a1;
          ndl_bond^[i].a2         := bond^[i].a2;
          ndl_bond^[i].btype      := bond^[i].btype;
          ndl_bond^[i].arom       := bond^[i].arom;
          ndl_bond^[i].ring_count := bond^[i].ring_count;  // new in v0.3d
          ndl_bond^[i].topo       := bond^[i].topo;        // new in v0.3d
          ndl_bond^[i].stereo     := bond^[i].stereo;      // new in v0.3d
        end;
    end;
  if (n_rings > 0) then
    begin
      for i := 1 to n_rings do
        begin
          for j := 1 to max_ringsize do ndl_ring^[i,j] := ring^[i,j];
        end;
      for i := 1 to max_rings do   // new in v0.3
        begin
          ndl_ringprop^[i].size     := ringprop^[i].size;
          ndl_ringprop^[i].arom     := ringprop^[i].arom;
          ndl_ringprop^[i].envelope := ringprop^[i].envelope;
        end;
    end;
  with ndl_molstat do
    begin
      n_QA := molstat.n_QA; n_QB := molstat.n_QB; n_chg := molstat.n_chg;
      n_C1 := molstat.n_C1; n_C2 := molstat.n_C2; n_C := molstat.n_C;
      n_CHB1p := molstat.n_CHB1p; n_CHB2p := molstat.n_CHB2p;
      n_CHB3p := molstat.n_CHB3p; n_CHB4 := molstat.n_CHB4;
      n_O2 := molstat.n_O2; n_O3  := molstat.n_O3;
      n_N1 := molstat.n_N1; n_N2 := molstat.n_N2; n_N3 := molstat.n_N3;
      n_S := molstat.n_S; n_SeTe := molstat.n_SeTe;
      n_F := molstat.n_F; n_Cl := molstat.n_Cl; n_Br := molstat.n_Br; n_I := molstat.n_I;
      n_P := molstat.n_P; n_B := molstat.n_B;
      n_Met := molstat.n_Met; n_X := molstat.n_X;
      n_b1 := molstat.n_b1; n_b2 := molstat.n_b2; n_b3 := molstat.n_b3; n_bar := molstat.n_bar;
      n_C1O := molstat.n_C1O; n_C2O := molstat.n_C2O; n_CN := molstat.n_CN; n_XY := molstat.n_XY;
      n_r3 := molstat.n_r3; n_r4 := molstat.n_r4; n_r5 := molstat.n_r5; n_r6 := molstat.n_r6;
      n_r7 := molstat.n_r7; n_r8 := molstat.n_r8; n_r9 := molstat.n_r9; n_r10 := molstat.n_r10;
      n_r11 := molstat.n_r11; n_r12 := molstat.n_r12; n_r13p := molstat.n_r13p;
      n_rN := molstat.n_rN; n_rN1 := molstat.n_rN1; n_rN2 := molstat.n_rN2; n_rN3p := molstat.n_rN3p;
      n_rO := molstat.n_rO; n_rO1 := molstat.n_rO1; n_rO2p := molstat.n_rO2p;
      n_rS := molstat.n_rS; n_rX := molstat.n_rX;
      n_rar := molstat.n_rar; n_rBz := molstat.n_rBz;  // v0.3l
      n_br2p := molstat.n_br2p;  // v0.3n
      {$IFDEF extended_molstat}                        // v0.3m
      n_psg01 := molstat.n_psg01; n_psg02 := molstat.n_psg02; n_psg13 := molstat.n_psg13; 
      n_psg14 := molstat.n_psg14; n_psg15 := molstat.n_psg15; n_psg16 := molstat.n_psg16; 
      n_psg17 := molstat.n_psg17; n_psg18 := molstat.n_psg18; 
      n_pstm := molstat.n_pstm; n_psla := molstat.n_psla; 
      {$ENDIF}
    end;
  // make sure some modes can be switched on only by the query file
  // and not by subsequent haystack file(s)
  if ez_flag   then ez_search := true;    // new in v0.3f
  if chir_flag then rs_search := true;    // new in v0.3f
  ndl_querymol := found_querymol;  // v0.3p
end;

// an alternative version (slightly slower with fpc 1.0.11 on Linux):
(*
procedure copy_mol_to_needle;
begin
  if (n_atoms = 0) then exit;
  try
    getmem(ndl_atom,n_atoms*sizeof(atom_rec));
    getmem(ndl_bond,n_bonds*sizeof(bond_rec));
    getmem(ndl_ring,sizeof(ringlist));
    getmem(ndl_ringprop,sizeof(ringprop_type));
  except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end;
  ndl_n_atoms := n_atoms;
  ndl_n_bonds := n_bonds;
  ndl_n_rings := n_rings;
  ndl_n_heavyatoms := n_heavyatoms;
  ndl_n_trueheavyatoms := n_trueheavyatoms;  // v0.4b
  ndl_n_heavybonds := n_heavybonds;
  ndl_molname := molname;
  ndl_n_Ctot := n_Ctot;
  ndl_n_Otot := n_Otot;
  ndl_n_Ntot := n_Ntot;
  move(atom^,ndl_atom^,n_atoms*sizeof(atom_rec));
  if (n_bonds > 0) then move(bond^,ndl_bond^,n_bonds*sizeof(bond_rec));
  if (n_rings > 0) then
    begin
      move(ring^,ndl_ring^,sizeof(ringlist));
      move(ringprop^,ndl_ringprop^,sizeof(ringprop_type));
    end;
  move(molstat,ndl_molstat,sizeof(molstat_rec));
  // make sure some modes can be switched on only by the query file
  // and not by subsequent haystack file(s)
  if ez_flag   then ez_search := true;    // new in v0.3f
  if chir_flag then rs_search := true;    // new in v0.3f
  ndl_querymol := found_querymol;         // v0.3p
end;
*)

procedure copy_mol_to_tmp;
var
  i, j : integer;
begin
  if (n_atoms = 0) then exit;
  try
    getmem(tmp_atom,n_atoms*sizeof(atom_rec));
    getmem(tmp_bond,n_bonds*sizeof(bond_rec));
    getmem(tmp_ring,sizeof(ringlist));
    getmem(tmp_ringprop,sizeof(ringprop_type));
  except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end;
  tmp_n_atoms := n_atoms;
  tmp_n_bonds := n_bonds;
  tmp_n_rings := n_rings;
  tmp_n_heavyatoms := n_heavyatoms;
  tmp_n_trueheavyatoms := n_trueheavyatoms;  // v0.4b
  tmp_n_heavybonds := n_heavybonds;
  tmp_molname := molname;
  tmp_n_Ctot := n_Ctot;
  tmp_n_Otot := n_Otot;
  tmp_n_Ntot := n_Ntot;
  for i := 1 to n_atoms do
    begin
      tmp_atom^[i].element        := atom^[i].element;
      tmp_atom^[i].atype          := atom^[i].atype;
      tmp_atom^[i].x              := atom^[i].x;
      tmp_atom^[i].y              := atom^[i].y;
      tmp_atom^[i].z              := atom^[i].z;
      tmp_atom^[i].formal_charge  := atom^[i].formal_charge;
      tmp_atom^[i].real_charge    := atom^[i].real_charge;
      tmp_atom^[i].Hexp           := atom^[i].Hexp;
      tmp_atom^[i].Htot           := atom^[i].Htot;
      tmp_atom^[i].neighbor_count := atom^[i].neighbor_count;
      tmp_atom^[i].ring_count     := atom^[i].ring_count;
      tmp_atom^[i].arom           := atom^[i].arom;
      tmp_atom^[i].stereo_care    := atom^[i].stereo_care;
      tmp_atom^[i].heavy          := atom^[i].heavy;  // v0.3l
      tmp_atom^[i].metal          := atom^[i].metal;  // v0.3l
      tmp_atom^[i].tag            := atom^[i].tag;  // v0.3o
      tmp_atom^[i].nucleon_number := atom^[i].nucleon_number;  // v0.3p
      tmp_atom^[i].radical_type   := atom^[i].radical_type;  // v0.3p
    end;
  if (n_bonds > 0) then
    begin
      for i := 1 to n_bonds do
        begin
          tmp_bond^[i].a1         := bond^[i].a1;
          tmp_bond^[i].a2         := bond^[i].a2;
          tmp_bond^[i].btype      := bond^[i].btype;
          tmp_bond^[i].arom       := bond^[i].arom;
          tmp_bond^[i].ring_count := bond^[i].ring_count;  // new in v0.3d
          tmp_bond^[i].topo       := bond^[i].topo;        // new in v0.3d
          tmp_bond^[i].stereo     := bond^[i].stereo;      // new in v0.3d
        end;
    end;
  if (n_rings > 0) then
    begin
      for i := 1 to n_rings do
        begin
          for j := 1 to max_ringsize do tmp_ring^[i,j] := ring^[i,j];
        end;
      for i := 1 to max_rings do   // new in v0.3
        begin
          tmp_ringprop^[i].size     := ringprop^[i].size;
          tmp_ringprop^[i].arom     := ringprop^[i].arom;
          tmp_ringprop^[i].envelope := ringprop^[i].envelope;
        end;
    end;
  with tmp_molstat do
    begin
      n_QA := molstat.n_QA; n_QB := molstat.n_QB; n_chg := molstat.n_chg;
      n_C1 := molstat.n_C1; n_C2 := molstat.n_C2; n_C := molstat.n_C;
      n_CHB1p := molstat.n_CHB1p; n_CHB2p := molstat.n_CHB2p;
      n_CHB3p := molstat.n_CHB3p; n_CHB4 := molstat.n_CHB4;
      n_O2 := molstat.n_O2; n_O3  := molstat.n_O3;
      n_N1 := molstat.n_N1; n_N2 := molstat.n_N2; n_N3 := molstat.n_N3;
      n_S := molstat.n_S; n_SeTe := molstat.n_SeTe;
      n_F := molstat.n_F; n_Cl := molstat.n_Cl; n_Br := molstat.n_Br; n_I := molstat.n_I;
      n_P := molstat.n_P; n_B := molstat.n_B;
      n_Met := molstat.n_Met; n_X := molstat.n_X;
      n_b1 := molstat.n_b1; n_b2 := molstat.n_b2; n_b3 := molstat.n_b3; n_bar := molstat.n_bar;
      n_C1O := molstat.n_C1O; n_C2O := molstat.n_C2O; n_CN := molstat.n_CN; n_XY := molstat.n_XY;
      n_r3 := molstat.n_r3; n_r4 := molstat.n_r4; n_r5 := molstat.n_r5; n_r6 := molstat.n_r6;
      n_r7 := molstat.n_r7; n_r8 := molstat.n_r8; n_r9 := molstat.n_r9; n_r10 := molstat.n_r10;
      n_r11 := molstat.n_r11; n_r12 := molstat.n_r12; n_r13p := molstat.n_r13p;
      n_rN := molstat.n_rN; n_rN1 := molstat.n_rN1; n_rN2 := molstat.n_rN2; n_rN3p := molstat.n_rN3p;
      n_rO := molstat.n_rO; n_rO1 := molstat.n_rO1; n_rO2p := molstat.n_rO2p;
      n_rS := molstat.n_rS; n_rX := molstat.n_rX;
      n_rar := molstat.n_rar; n_rBz := molstat.n_rBz;  // v0.3l
      n_br2p := molstat.n_br2p;  // v0.3n
      {$IFDEF extended_molstat}                        // v0.3m
      n_psg01 := molstat.n_psg01; n_psg02 := molstat.n_psg02; n_psg13 := molstat.n_psg13; 
      n_psg14 := molstat.n_psg14; n_psg15 := molstat.n_psg15; n_psg16 := molstat.n_psg16; 
      n_psg17 := molstat.n_psg17; n_psg18 := molstat.n_psg18; 
      n_pstm := molstat.n_pstm; n_psla := molstat.n_psla; 
      {$ENDIF}
    end;
  // make sure some modes can be switched on only by the query file
  // and not by subsequent haystack file(s)
  if ez_flag   then ez_search := true;    // new in v0.3f
  if chir_flag then rs_search := true;    // new in v0.3f
end;

// an alternative version (slightly slower with fpc 1.0.11 on Linux):
(*
procedure copy_mol_to_tmp;
begin
  if (n_atoms = 0) then exit;
  try
    getmem(tmp_atom,n_atoms*sizeof(atom_rec));
    getmem(tmp_bond,n_bonds*sizeof(bond_rec));
    getmem(tmp_ring,sizeof(ringlist));
    getmem(tmp_ringprop,sizeof(ringprop_type));
  except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end;
  tmp_n_atoms := n_atoms;
  tmp_n_bonds := n_bonds;
  tmp_n_rings := n_rings;
  tmp_n_heavyatoms := n_heavyatoms;
  tmp_n_trueheavyatoms := n_trueheavyatoms;  // v0.4b
  tmp_n_heavybonds := n_heavybonds;
  tmp_molname := molname;
  tmp_n_Ctot := n_Ctot;
  tmp_n_Otot := n_Otot;
  tmp_n_Ntot := n_Ntot;
  move(atom^,tmp_atom^,n_atoms*sizeof(atom_rec));
  if (n_bonds > 0) then move(bond^,tmp_bond^,n_bonds*sizeof(bond_rec));
  if (n_rings > 0) then
    begin
      move(ring^,tmp_ring^,sizeof(ringlist));
      move(ringprop^,tmp_ringprop^,sizeof(ringprop_type));
    end;
  move(molstat,tmp_molstat,sizeof(molstat_rec));
  // make sure some modes can be switched on only by the query file
  // and not by subsequent haystack file(s)
  if ez_flag   then ez_search := true;    // new in v0.3f
  if chir_flag then rs_search := true;    // new in v0.3f
end;
*)

procedure copy_tmp_to_mol;
var
  i, j : integer;
begin
  if (tmp_n_atoms = 0) then exit;
  n_atoms := tmp_n_atoms;
  n_bonds := tmp_n_bonds;
  n_rings := tmp_n_rings;
  n_heavyatoms := tmp_n_heavyatoms;
  n_trueheavyatoms := tmp_n_trueheavyatoms;  // v0.4b
  n_heavybonds := tmp_n_heavybonds;
  molname := tmp_molname;
  n_Ctot := tmp_n_Ctot;
  n_Otot := tmp_n_Otot;
  n_Ntot := tmp_n_Ntot;
  try
    getmem(atom,n_atoms*sizeof(atom_rec));
    getmem(bond,n_bonds*sizeof(bond_rec));
    getmem(ring,sizeof(ringlist));
    getmem(ringprop,sizeof(ringprop_type));
  except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end;
  for i := 1 to tmp_n_atoms do
    begin
      atom^[i].element        := tmp_atom^[i].element;
      atom^[i].atype          := tmp_atom^[i].atype;
      atom^[i].x              := tmp_atom^[i].x;
      atom^[i].y              := tmp_atom^[i].y;
      atom^[i].z              := tmp_atom^[i].z;
      atom^[i].formal_charge  := tmp_atom^[i].formal_charge;
      atom^[i].real_charge    := tmp_atom^[i].real_charge;
      atom^[i].Hexp           := tmp_atom^[i].Hexp;
      atom^[i].Htot           := tmp_atom^[i].Htot;
      atom^[i].neighbor_count := tmp_atom^[i].neighbor_count;
      atom^[i].ring_count     := tmp_atom^[i].ring_count;
      atom^[i].arom           := tmp_atom^[i].arom;
      atom^[i].stereo_care    := tmp_atom^[i].stereo_care;
      atom^[i].heavy          := tmp_atom^[i].heavy;  // v0.3l
      atom^[i].metal          := tmp_atom^[i].metal;  // v0.3l
      atom^[i].tag            := tmp_atom^[i].tag;  // v0.3o
      atom^[i].nucleon_number := tmp_atom^[i].nucleon_number;  // v0.3p
      atom^[i].radical_type   := tmp_atom^[i].radical_type;  // v0.3p
    end;
  if (tmp_n_bonds > 0) then
    begin
      for i := 1 to tmp_n_bonds do
        begin
          bond^[i].a1         := tmp_bond^[i].a1;
          bond^[i].a2         := tmp_bond^[i].a2;
          bond^[i].btype      := tmp_bond^[i].btype;
          bond^[i].arom       := tmp_bond^[i].arom;
          bond^[i].ring_count := tmp_bond^[i].ring_count;  // new in v0.3d
          bond^[i].topo       := tmp_bond^[i].topo;        // new in v0.3d
          bond^[i].stereo     := tmp_bond^[i].stereo;      // new in v0.3d
        end;
    end;
  if (tmp_n_rings > 0) then
    begin
      for i := 1 to tmp_n_rings do
        begin
          for j := 1 to max_ringsize do ring^[i,j] := tmp_ring^[i,j];
        end;
      for i := 1 to max_rings do   // new in v0.3
        begin
          ringprop^[i].size     := tmp_ringprop^[i].size;
          ringprop^[i].arom     := tmp_ringprop^[i].arom;
          ringprop^[i].envelope := tmp_ringprop^[i].envelope;
        end;
    end;
  with molstat do
    begin
      n_QA := tmp_molstat.n_QA; n_QB := tmp_molstat.n_QB; n_chg := tmp_molstat.n_chg;
      n_C1 := tmp_molstat.n_C1; n_C2 := tmp_molstat.n_C2; n_C := tmp_molstat.n_C;
      n_CHB1p := tmp_molstat.n_CHB1p; n_CHB2p := tmp_molstat.n_CHB2p;
      n_CHB3p := tmp_molstat.n_CHB3p; n_CHB4 := tmp_molstat.n_CHB4;
      n_O2 := tmp_molstat.n_O2; n_O3  := tmp_molstat.n_O3;
      n_N1 := tmp_molstat.n_N1; n_N2 := tmp_molstat.n_N2; n_N3 := tmp_molstat.n_N3;
      n_S := tmp_molstat.n_S; n_SeTe := tmp_molstat.n_SeTe;
      n_F := tmp_molstat.n_F; n_Cl := tmp_molstat.n_Cl; n_Br := tmp_molstat.n_Br; n_I := tmp_molstat.n_I;
      n_P := tmp_molstat.n_P; n_B := tmp_molstat.n_B;
      n_Met := tmp_molstat.n_Met; n_X := tmp_molstat.n_X;
      n_b1 := tmp_molstat.n_b1; n_b2 := tmp_molstat.n_b2; n_b3 := tmp_molstat.n_b3; n_bar := tmp_molstat.n_bar;
      n_C1O := tmp_molstat.n_C1O; n_C2O := tmp_molstat.n_C2O; n_CN := tmp_molstat.n_CN; n_XY := tmp_molstat.n_XY;
      n_r3 := tmp_molstat.n_r3; n_r4 := tmp_molstat.n_r4; n_r5 := tmp_molstat.n_r5; n_r6 := tmp_molstat.n_r6;
      n_r7 := tmp_molstat.n_r7; n_r8 := tmp_molstat.n_r8; n_r9 := tmp_molstat.n_r9; n_r10 := tmp_molstat.n_r10;
      n_r11 := tmp_molstat.n_r11; n_r12 := tmp_molstat.n_r12; n_r13p := tmp_molstat.n_r13p;
      n_rN := tmp_molstat.n_rN; n_rN1 := tmp_molstat.n_rN1; n_rN2 := tmp_molstat.n_rN2; n_rN3p := tmp_molstat.n_rN3p;
      n_rO := tmp_molstat.n_rO; n_rO1 := tmp_molstat.n_rO1; n_rO2p := tmp_molstat.n_rO2p;
      n_rS := tmp_molstat.n_rS; n_rX := tmp_molstat.n_rX;
      n_rar := tmp_molstat.n_rar; n_rBz := tmp_molstat.n_rBz;  // v0.3l
      n_br2p := tmp_molstat.n_br2p;  // v0.3n
      {$IFDEF extended_molstat}
      n_psg01 := tmp_molstat.n_psg01; n_psg02 := tmp_molstat.n_psg02; n_psg13 := tmp_molstat.n_psg13; 
      n_psg14 := tmp_molstat.n_psg14; n_psg15 := tmp_molstat.n_psg15; n_psg16 := tmp_molstat.n_psg16; 
      n_psg17 := tmp_molstat.n_psg17; n_psg18 := tmp_molstat.n_psg18; 
      n_pstm := tmp_molstat.n_pstm; n_psla := tmp_molstat.n_psla; 
      {$ENDIF}
    end;
  // make sure some modes can be switched on only by the query file
  // and not by subsequent haystack file(s)
  if ez_flag   then ez_search := true;
  if chir_flag then rs_search := true;
end;

// an alternative version (slightly slower with fpc 1.0.11 on Linux):
(*
procedure copy_tmp_to_mol;
begin
  if (n_atoms = 0) then exit;
  n_atoms := tmp_n_atoms;
  n_bonds := tmp_n_bonds;
  n_rings := tmp_n_rings;
  n_heavyatoms := tmp_n_heavyatoms;
  n_trueheavyatoms := tmp_n_trueheavyatoms;  // v0.4b
  n_heavybonds := tmp_n_heavybonds;
  molname := tmp_molname;
  n_Ctot := tmp_n_Ctot;
  n_Otot := tmp_n_Otot;
  n_Ntot := tmp_n_Ntot;
  try
    getmem(atom,n_atoms*sizeof(atom_rec));
    getmem(bond,n_bonds*sizeof(bond_rec));
    getmem(ring,sizeof(ringlist));
    getmem(ringprop,sizeof(ringprop_type));
  except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end;
  move(tmp_atom^,atom^,tmp_n_atoms*sizeof(atom_rec));
  if (tmp_n_bonds > 0) then move(tmp_bond^,bond^,tmp_n_bonds*sizeof(bond_rec));
  if (tmp_n_rings > 0) then
    begin
      move(tmp_ring^,ring^,sizeof(ringlist));
      move(tmp_ringprop^,ringprop^,sizeof(ringprop_type));
    end;
  move(tmp_molstat,molstat,sizeof(molstat_rec));
  // make sure some modes can be switched on only by the query file
  // and not by subsequent haystack file(s)
  if ez_flag   then ez_search := true;
  if chir_flag then rs_search := true;
end;
*)


procedure get_ringstat(r_id:integer);
var
  i,j  : integer;
  testring : ringpath_type;
  ring_size : integer;
  a_ref : integer;
  elem : str2;
  nN, nO, nS, nX : integer;
begin
  nN := 0; nO := 0; nS := 0; nX := 0;
  if (r_id < 1) or (r_id > n_rings) then exit;
  fillchar(testring,sizeof(ringpath_type),0);
  ring_size := ringprop^[r_id].size;  // v0.3j
  for j := 1 to ring_size do testring[j] := ring^[r_id,j];  // v0.3j
  {$IFDEF reduced_SAR}  
  if (ring_size > 2) and (ringprop^[r_id].envelope = false) then  // v0.3n: ignore envelope rings
  {$ELSE}
  if (ring_size > 2) then
  {$ENDIF}
    begin
      for i := 1 to ring_size do
        begin
          a_ref := testring[i];
          elem := atom^[a_ref].element;
          if (elem <> 'C ') and (elem <> 'A ') then 
            begin
              inc(nX);   // general heteroatom count
              if elem = 'N ' then inc(nN);
              if elem = 'O ' then inc(nO);
              if elem = 'S ' then inc(nS);
            end;
        end;
      if nN > 0 then
        begin
          inc(molstat.n_rN);
          if nN = 1 then inc(molstat.n_rN1);
          if nN = 2 then inc(molstat.n_rN2);
          if nN > 2 then inc(molstat.n_rN3p);
        end;  
      if nO > 0 then
        begin
          inc(molstat.n_rO);
          if nO = 1 then inc(molstat.n_rO1);
          if nO = 2 then inc(molstat.n_rO2p);
        end;
      if nS > 0 then inc(molstat.n_rS);
      if nX > 0 then inc(molstat.n_rX);
      // general ringsize descriptors; v0.3m
      case ring_size of
        3 : inc(molstat.n_r3);
        4 : inc(molstat.n_r4);
        5 : inc(molstat.n_r5);
        6 : inc(molstat.n_r6);
        7 : inc(molstat.n_r7);
        8 : inc(molstat.n_r8);
        9 : inc(molstat.n_r9);
        10 : inc(molstat.n_r10);
        11 : inc(molstat.n_r11);
        12 : inc(molstat.n_r12);
        else inc(molstat.n_r13p);
      end;  // end v0.3m       
    end;
end;


procedure get_molstat;
var
  i : integer;
  elem  : str2;
  atype : str3;
  a1, a2 : integer;
  a1el, a2el : str2;
  btype : char;
  hbc : integer;
  n_b2formal : integer;  // new in v0.2e
begin
  if n_atoms = 0 then exit;
  with molstat do
    begin
      for i := 1 to n_atoms do
        begin
          if atom^[i].heavy then
            begin
              elem  := atom^[i].element;
              atype := atom^[i].atype;
              if (atype = 'C1 ') then inc(n_C1);
              if (atype = 'C2 ') or (atype = 'CAR') then inc(n_C2);
              if (elem  = 'C ') then inc(n_C);
              if (atype = 'O2 ') then inc(n_O2);
              if (atype = 'O3 ') then inc(n_O3);
              if (atype = 'N1 ') then inc(n_N1);
              if (atype = 'N2 ') or (atype = 'NAR') or 
                ((atype = 'NAM') and (atom^[i].arom = true)) then inc(n_N2);  // v0.3n
              if (atype = 'N3 ') or (atype = 'NPL') or (atype = 'N3+') or 
                ((atype = 'NAM') and (atom^[i].arom = false)) then inc(n_N3);  // v0.3n
              if (elem = 'A ') then inc(n_QA);  // query atom
              if (elem = 'Q ') then inc(n_QA);  // query atom
              if (elem = 'X ') then inc(n_QA);  // query atom (halogen)  // v0.3p
              if (elem = 'S ') then inc(n_S);
              if (elem = 'SE') then inc(n_SeTe);
              if (elem = 'TE') then inc(n_SeTe);
              if (elem = 'F ') then inc(n_F);
              if (elem = 'CL') then inc(n_Cl);
              if (elem = 'BR') then inc(n_Br);
              if (elem = 'I ') then inc(n_I);
              if (elem = 'P ') then inc(n_P);
              if (elem = 'B ') then inc(n_B);
	      // check for known metals
	      if (atom^[i].metal) then inc(n_Met);  // v0.3l
              // still missing: unknown elements

              // check number of heteroatom bonds per C atom
              if elem = 'C ' then
                begin
                  hbc := raw_hetbond_count(i);   // new in v0.2j (replaces hetbond_count)
                  if hbc >= 1 then inc(n_CHB1p);
                  if hbc >= 2 then inc(n_CHB2p);
                  if hbc >= 3 then inc(n_CHB3p);
                  if hbc  = 4 then inc(n_CHB4);
                end;
              if atom^[i].formal_charge <> 0 then 
                begin
                  inc(n_chg);
                  inc(n_charges);
                end;
	      // check for "other" elements;  v0.3l
              if not (atom^[i].metal or (elem = 'C ') or (elem = 'N ') or (elem = 'O ') or
                //(elem = 'F ') or (elem = 'CL') or (elem = 'BR') or (elem = 'I ') or  // leave halogens as type X, v0.3m
                (elem = 'S ') or (elem = 'SE') or (elem = 'TE') or (elem = 'P ') or
                (elem = 'B ') or (elem = 'A ') or (elem = 'Q ') ) then inc(n_X);
              {$IFDEF extended_molstat}
              if (elem = 'LI') or (elem = 'NA') or (elem = 'K ') or (elem = 'RB') or
                 (elem = 'CS') or (elem = 'FR') then inc(n_psg01); 
              if (elem = 'BE') or (elem = 'MG') or (elem = 'CA') or (elem = 'SR') or
                 (elem = 'BA') or (elem = 'RA') then inc(n_psg02); 
              if (elem = 'B ') or (elem = 'AL') or (elem = 'GA') or (elem = 'IN') or
                 (elem = 'TL') then inc(n_psg13); 
              if (elem = 'C ') or (elem = 'SI') or (elem = 'GE') or (elem = 'SN') or
                 (elem = 'PB') then inc(n_psg14); 
              if (elem = 'N ') or (elem = 'P ') or (elem = 'AS') or (elem = 'SB') or
                 (elem = 'BI') then inc(n_psg15); 
              if (elem = 'O ') or (elem = 'S ') or (elem = 'SE') or (elem = 'TE') or
                 (elem = 'PO') then inc(n_psg16); 
              if (elem = 'F ') or (elem = 'CL') or (elem = 'BR') or (elem = 'I ') or
                 (elem = 'AT') then inc(n_psg17); 
              if (elem = 'HE') or (elem = 'NE') or (elem = 'AR') or (elem = 'KR') or
                 (elem = 'XE') or (elem = 'RN') then inc(n_psg18); 
              if (elem = 'SC') or (elem = 'Y ') or (elem = 'LU') or (elem = 'LR') or
                 (elem = 'TI') or (elem = 'ZR') or (elem = 'HF') or (elem = 'RF') or
                 (elem = 'V ') or (elem = 'NB') or (elem = 'TA') or (elem = 'DB') or
                 (elem = 'CR') or (elem = 'MO') or (elem = 'W ') or (elem = 'SG') or
                 (elem = 'MN') or (elem = 'TC') or (elem = 'RE') or (elem = 'BH') or
                 (elem = 'FE') or (elem = 'RU') or (elem = 'OS') or (elem = 'HS') or
                 (elem = 'CO') or (elem = 'RH') or (elem = 'IR') or (elem = 'MT') or
                 (elem = 'NI') or (elem = 'PD') or (elem = 'PT') or (elem = 'DS') or
                 (elem = 'CU') or (elem = 'AG') or (elem = 'AU') or (elem = 'RG') or
                 (elem = 'ZN') or (elem = 'CD') or (elem = 'HG') then inc(n_pstm); 
              if (elem = 'LA') or (elem = 'CE') or (elem = 'PR') or (elem = 'ND') or
                 (elem = 'PM') or (elem = 'SM') or (elem = 'EU') or (elem = 'GD') or
                 (elem = 'TB') or (elem = 'DY') or (elem = 'HO') or (elem = 'ER') or
                 (elem = 'TM') or (elem = 'YB') or
                 (elem = 'AC') or (elem = 'TH') or (elem = 'PA') or (elem = 'U ') or
                 (elem = 'NP') or (elem = 'PU') or (elem = 'AM') or (elem = 'CM') or
                 (elem = 'BK') or (elem = 'CF') or (elem = 'ES') or (elem = 'FM') or
                 (elem = 'MD') or (elem = 'NO') then inc(n_psla); 
              {$ENDIF}
            end;  // is heavy
        end;  // atoms
      if n_bonds > 0 then
        begin
          for i := 1 to n_bonds do
            begin
              a1 := bond^[i].a1; a2 := bond^[i].a2;
              a1el := atom^[a1].element;
              a2el := atom^[a2].element;
              btype := bond^[i].btype;
              if bond^[i].arom then inc(n_bar) else
                begin                  // v0.3n: ignore bonds to (explicit) hydrogens
                  if (btype = 'S') and (atom^[a1].heavy and atom^[a2].heavy) then inc(n_b1);
                  if (btype = 'D') then inc(n_b2);
                  if (btype = 'T') then inc(n_b3);
                end;
              if ((a1el = 'C ') and (a2el = 'O ')) or ((a1el = 'O ') and (a2el = 'C ')) then
                begin
                  if (btype = 'S') then inc(n_C1O);
                  if (btype = 'D') then inc(n_C2O);
                end;
              if ((a1el = 'C ') and (a2el = 'N ')) or ((a1el = 'N ') and (a2el = 'C ')) then inc(n_CN);
              if ((a1el <> 'C ') and (atom^[a1].heavy) and ((a2el <> 'C ') and (atom^[a2].heavy))) then inc(n_XY);
              // new in v0.3n: number of bonds belonging to more than one ring
              if (bond^[i].ring_count > 1) then inc(n_br2p);
            end;
        end; // bonds
      if n_rings > 0 then
        begin
          n_b2formal := 0;        // v0.3n
          n_countablerings := 0;  // v0.3n
          for i := 1 to n_rings do
            begin
              if (ringprop^[i].envelope = false) then inc(n_countablerings);  // v0.3n
              if (is_arene(i) and (ringprop^[i].envelope = false)) then   // v0.3n: ignore envelope rings
                begin
                  inc(n_rar);
                  if ((ringprop^[i].size = 6) and (is_heterocycle(i) = false)) then inc(n_rBz);  // v0.3l
                end;
              get_ringstat(i);
              if ((ringprop^[i].arom = true) and (ringprop^[i].envelope = false)) then 
                inc(n_b2formal);  // new in v0.3n; replaces assignment below
            end;
          //n_b2formal := n_rar;  // new in v0.2e; adds 1 formal double bond for each aromatic ring
                                // in order to allow an isolated double bond in the needle
                                // to be matched as a ring fragment of an aromatic ring
          if (n_b2formal > (n_bar div 2)) then n_b2formal := n_bar div 2;
          n_b2 := n_b2 + n_b2formal;
        end; // rings
    end;
end;

procedure fix_ssr_ringcounts;  // new in v0.3n
  // if SAR -> SSR fallback happens, set some molstat values
  // to a maximum (ring counts for various ring sizes);
  // this should be necessary only for ring sizes which
  // are a) too large for the SSR (depending on ssr_vringsize)
  // and b) which are likely to contain "envelope rings"
  // (size 6 and above)
begin
//  if (molstat.n_r3 = 0) then molstat.n_r3 := max_rings;
//  if (molstat.n_r4 = 0) then molstat.n_r4 := max_rings;
//  if (molstat.n_r5 = 0) then molstat.n_r5 := max_rings;
  if (molstat.n_r6 = 0) then molstat.n_r6 := max_rings;
  if (molstat.n_r7 = 0) then molstat.n_r7 := max_rings;
  if (molstat.n_r8 = 0) then molstat.n_r8 := max_rings;
  if (molstat.n_r9 = 0) then molstat.n_r9 := max_rings;
  if (molstat.n_r10 = 0) then molstat.n_r10 := max_rings;
  if (molstat.n_r11 = 0) then molstat.n_r11 := max_rings;
  if (molstat.n_r12 = 0) then molstat.n_r12 := max_rings;
  if (molstat.n_r13p =0) then molstat.n_r13p := max_rings;
end;



procedure write_molstat;
begin
  if auto_ssr then fix_ssr_ringcounts;  // v0.3n
  with molstat do
    begin
      write('n_atoms:',n_heavyatoms,';');  // count only non-H atoms (some molfiles contain explicit H's)
      if (opt_molstat_v = false) then  // v0.4d
        begin
          if n_bonds > 0 then write('n_bonds:',n_heavybonds,';');  // count only bonds between non-H atoms
          {$IFDEF reduced_SAR}
          if n_rings > 0 then write('n_rings:',n_countablerings,';');  // changed to non-envelope rings in v0.3n
          {$ELSE}
          if n_rings > 0 then write('n_rings:',n_rings,';');  // changed to non-envelope rings in v0.3n
          {$ENDIF}
          //      if n_QA    > 0 then write('n_QA:',n_QA,';');
          //      if n_QB    > 0 then write('n_QB:',n_QB,';');
          if opt_chg and (n_chg   > 0) then write('n_chg:',n_chg,';');   // v0.3p
          if n_C1    > 0 then write('n_C1:',n_C1,';');
          if n_C2    > 0 then write('n_C2:',n_C2,';');
          // requirement of a given number of sp3 carbons might be too restrictive,
          // so we use the total number of carbons instead  (initially used variable n_C3 is now n_C)
          if n_C    > 0 then write('n_C:',n_C,';');
          if n_CHB1p > 0 then write('n_CHB1p:',n_CHB1p,';');
          if n_CHB2p > 0 then write('n_CHB2p:',n_CHB2p,';');
          if n_CHB3p > 0 then write('n_CHB3p:',n_CHB3p,';');
          if n_CHB4  > 0 then write('n_CHB4:',n_CHB4,';');
          if n_O2    > 0 then write('n_O2:',n_O2,';');
          if n_O3    > 0 then write('n_O3:',n_O3,';');
          if n_N1    > 0 then write('n_N1:',n_N1,';');
          if n_N2    > 0 then write('n_N2:',n_N2,';');
          if n_N3    > 0 then write('n_N3:',n_N3,';');
          if n_S     > 0 then write('n_S:',n_S,';');
          if n_SeTe  > 0 then write('n_SeTe:',n_SeTe,';');
          if n_F     > 0 then write('n_F:',n_F,';');
          if n_Cl    > 0 then write('n_Cl:',n_Cl,';');
          if n_Br    > 0 then write('n_Br:',n_Br,';');
          if n_I     > 0 then write('n_I:',n_I,';');
          if n_P     > 0 then write('n_P:',n_P,';');
          if n_B     > 0 then write('n_B:',n_B,';');
          if n_Met   > 0 then write('n_Met:',n_Met,';');
          if n_X     > 0 then write('n_X:',n_X,';');
          if n_b1    > 0 then write('n_b1:',n_b1,';');
          if n_b2    > 0 then write('n_b2:',n_b2,';');
          if n_b3    > 0 then write('n_b3:',n_b3,';');
          if n_bar   > 0 then write('n_bar:',n_bar,';');
          if n_C1O   > 0 then write('n_C1O:',n_C1O,';');
          if n_C2O   > 0 then write('n_C2O:',n_C2O,';');
          if n_CN    > 0 then write('n_CN:',n_CN,';');
          if n_XY    > 0 then write('n_XY:',n_XY,';');
          if n_r3    > 0 then write('n_r3:',n_r3,';');
          if n_r4    > 0 then write('n_r4:',n_r4,';');
          if n_r5    > 0 then write('n_r5:',n_r5,';');
          if n_r6    > 0 then write('n_r6:',n_r6,';');
          if n_r7    > 0 then write('n_r7:',n_r7,';');
          if n_r8    > 0 then write('n_r8:',n_r8,';');
          if n_r9    > 0 then write('n_r9:',n_r9,';');
          if n_r10   > 0 then write('n_r10:',n_r10,';');
          if n_r11   > 0 then write('n_r11:',n_r11,';');
          if n_r12   > 0 then write('n_r12:',n_r12,';');
          if n_r13p  > 0 then write('n_r13p:',n_r13p,';');
          if n_rN    > 0 then write('n_rN:',n_rN,';');
          if n_rN1   > 0 then write('n_rN1:',n_rN1,';');
          if n_rN2   > 0 then write('n_rN2:',n_rN2,';');
          if n_rN3p  > 0 then write('n_rN3p:',n_rN3p,';');
          if n_rO    > 0 then write('n_rO:',n_rO,';');
          if n_rO1   > 0 then write('n_rO1:',n_rO1,';');
          if n_rO2p  > 0 then write('n_rO2p:',n_rO2p,';');
          if n_rS    > 0 then write('n_rS:',n_rS,';');
          if n_rX    > 0 then write('n_rX:',n_rX,';');
          if n_rar   > 0 then write('n_rar:',n_rar,';');
          {$IFDEF extended_molstat}
          if n_rBz   > 0 then write('n_rbz:',n_rBz,';');
          if n_br2p  > 0 then write('n_br2p:',n_br2p,';');
          if n_psg01 > 0 then write('n_psg01:',n_psg01,';');
          if n_psg02 > 0 then write('n_psg02:',n_psg02,';');
          if n_psg13 > 0 then write('n_psg13:',n_psg13,';');
          if n_psg14 > 0 then write('n_psg14:',n_psg14,';');
          if n_psg15 > 0 then write('n_psg15:',n_psg15,';');
          if n_psg16 > 0 then write('n_psg16:',n_psg16,';');
          if n_psg17 > 0 then write('n_psg17:',n_psg17,';');
          if n_psg18 > 0 then write('n_psg18:',n_psg18,';');
          if n_pstm > 0 then write('n_pstm:',n_pstm,';');
          if n_psla > 0 then write('n_psla:',n_psla,';');
          {$ENDIF}
        end else
        begin   // full output, even if field is 0;  v0.4d
          write('n_bonds:',n_heavybonds,';');  // count only bonds between non-H atoms
          {$IFDEF reduced_SAR}
          write('n_rings:',n_countablerings,';');  // changed to non-envelope rings in v0.3n
          {$ELSE}
          write('n_rings:',n_rings,';');  // changed to non-envelope rings in v0.3n
          {$ENDIF}
          write('n_QA:',n_QA,';');
          write('n_QB:',n_QB,';');
          if opt_chg then write('n_chg:',n_chg,';');   // v0.3p
          write('n_C1:',n_C1,';');
          write('n_C2:',n_C2,';');
          write('n_C:',n_C,';');
          write('n_CHB1p:',n_CHB1p,';');
          write('n_CHB2p:',n_CHB2p,';');
          write('n_CHB3p:',n_CHB3p,';');
          write('n_CHB4:',n_CHB4,';');
          write('n_O2:',n_O2,';');
          write('n_O3:',n_O3,';');
          write('n_N1:',n_N1,';');
          write('n_N2:',n_N2,';');
          write('n_N3:',n_N3,';');
          write('n_S:',n_S,';');
          write('n_SeTe:',n_SeTe,';');
          write('n_F:',n_F,';');
          write('n_Cl:',n_Cl,';');
          write('n_Br:',n_Br,';');
          write('n_I:',n_I,';');
          write('n_P:',n_P,';');
          write('n_B:',n_B,';');
          write('n_Met:',n_Met,';');
          write('n_X:',n_X,';');
          write('n_b1:',n_b1,';');
          write('n_b2:',n_b2,';');
          write('n_b3:',n_b3,';');
          write('n_bar:',n_bar,';');
          write('n_C1O:',n_C1O,';');
          write('n_C2O:',n_C2O,';');
          write('n_CN:',n_CN,';');
          write('n_XY:',n_XY,';');
          write('n_r3:',n_r3,';');
          write('n_r4:',n_r4,';');
          write('n_r5:',n_r5,';');
          write('n_r6:',n_r6,';');
          write('n_r7:',n_r7,';');
          write('n_r8:',n_r8,';');
          write('n_r9:',n_r9,';');
          write('n_r10:',n_r10,';');
          write('n_r11:',n_r11,';');
          write('n_r12:',n_r12,';');
          write('n_r13p:',n_r13p,';');
          write('n_rN:',n_rN,';');
          write('n_rN1:',n_rN1,';');
          write('n_rN2:',n_rN2,';');
          write('n_rN3p:',n_rN3p,';');
          write('n_rO:',n_rO,';');
          write('n_rO1:',n_rO1,';');
          write('n_rO2p:',n_rO2p,';');
          write('n_rS:',n_rS,';');
          write('n_rX:',n_rX,';');
          write('n_rar:',n_rar,';');
          {$IFDEF extended_molstat}
          write('n_rbz:',n_rBz,';');
          write('n_br2p:',n_br2p,';');
          write('n_psg01:',n_psg01,';');
          write('n_psg02:',n_psg02,';');
          write('n_psg13:',n_psg13,';');
          write('n_psg14:',n_psg14,';');
          write('n_psg15:',n_psg15,';');
          write('n_psg16:',n_psg16,';');
          write('n_psg17:',n_psg17,';');
          write('n_psg18:',n_psg18,';');
          write('n_pstm:',n_pstm,';');
          write('n_psla:',n_psla,';');
          {$ENDIF}
        end;
      writeln;
    end;
end;


procedure write_molstat_X;
begin
  if auto_ssr then fix_ssr_ringcounts;  // v0.3n
  with molstat do
    begin
      {$IFDEF reduced_SAR}
      write(n_heavyatoms,','); write(n_heavybonds,','); write(n_countablerings,',');  // v0.3n: n_rings =?> n_countablerings
      {$ELSE}
      write(n_heavyatoms,','); write(n_heavybonds,','); write(n_rings,',');  // v0.3n: n_rings ==> n_countablerings
      {$ENDIF}
      write(n_QA,','); write(n_QB,','); 
      if opt_chg then write(n_chg,',') else write('0,');  // v0.3p
      write(n_C1,','); write(n_C2,','); write(n_C,',');
      write(n_CHB1p,','); write(n_CHB2p,','); write(n_CHB3p,','); write(n_CHB4,',');
      write(n_O2,','); write(n_O3,',');
      write(n_N1,','); write(n_N2,','); write(n_N3,',');
      write(n_S,','); write(n_SeTe,',');
      write(n_F,','); write(n_Cl,','); write(n_Br,','); write(n_I,',');
      write(n_P,','); write(n_B,',');
      write(n_Met,','); write(n_X,',');
      write(n_b1,','); write(n_b2,','); write(n_b3,','); write(n_bar,',');
      write(n_C1O,','); write(n_C2O,','); write(n_CN,','); write(n_XY,',');
      write(n_r3,','); write(n_r4,','); write(n_r5,','); write(n_r6,',');
      write(n_r7,','); write(n_r8,','); write(n_r9,','); write(n_r10,','); 
      write(n_r11,','); write(n_r12,','); write(n_r13p,',');
      write(n_rN,','); write(n_rN1,','); write(n_rN2,','); write(n_rN3p,',');
      write(n_rO,','); write(n_rO1,','); write(n_rO2p,',');
      write(n_rS,','); write(n_rX,',');
      write(n_rar);
      {$IFDEF extended_molstat}
      write(',',n_rBz);
      write(',',n_br2p);
      write(',',n_psg01);
      write(',',n_psg02);
      write(',',n_psg13);
      write(',',n_psg14);
      write(',',n_psg15);
      write(',',n_psg16);
      write(',',n_psg17);
      write(',',n_psg18);
      write(',',n_pstm);
      write(',',n_psla);
      {$ENDIF}
      writeln
    end;
end;


// routines for substructure matching


function find_ndl_ref_atom: integer;
var
  i : integer;
  score : integer;
  score_max : integer; // v0.4b
  index : integer;
  n_nb, n_hc : integer;
begin
  // finds a characteristic atom in the needle molecule,
  // i.e., one with as many substituents as possible and
  // with as many heteroatom substitutents as possible;
  // added in v0.2d: make sure that reference atom is a heavy atom
  // and not (accidentally) an explicit hydrogen;
  // new in v0.3d: special treatment in case of E/Z geometry search
  // to ensure that the entire A-B=C-D fragment is enclosed in one
  // matchpath, regardless where the recursive search starts;
  // refined in v0.3f: exclude only alkene-C as reference atoms
  // added in v0.3o: needle atom must be "tagged" in order to be
  // selected (prevents unconnected fragments from being overlooked)
  if ndl_n_atoms = 0 then exit;
  score := -1;
  score_max := -1;  // v0.4b
  index := 0;
  if (ez_search) and (ndl_n_heavyatoms > 2) then
    begin
      for i := 1 to ndl_n_atoms do
        begin  // ignore sp2-carbons if not aromatic
          //if ((ndl_atom^[i].atype <> 'C2 ') or (ndl_atom^[i].arom = true)) then
          if (ndl_alkene_C(i) = false) and (ndl_atom^[i].tag) then  // v0.3o
            begin
              n_nb := ndl_atom^[i].neighbor_count;
              n_hc := ndl_hetatom_count(i);
              if (11*n_nb + 7*n_hc > score) and (ndl_atom^[i].heavy) then  // v0.3j
                begin
                  index := i;
                  score := 11*n_nb + 7*n_hc;  // changed in v0.3j
                end;
            end;
          if (score > score_max) then score_max := score;  // v0.4b
        end;
    end;
  // it is possible that no suitable reference atom has been found here
  // (e.g., with "pure" polyenes), so we need a fallback option anyway
  if (index = 0) then
    begin
      ez_search := false;  // just in case it was true
      opt_geom := false;   // just in case it was true
      for i := 1 to ndl_n_atoms do
        begin
          n_nb := ndl_atom^[i].neighbor_count;
          n_hc := ndl_hetatom_count(i);
          if (11*n_nb + 7*n_hc > score) and (ndl_atom^[i].heavy) and  // v0.3j
             (ndl_atom^[i].tag) then  // v0.3o
            begin
              index := i;
              score := 11*n_nb + 7*n_hc;  // changed in v0.3j
            end;
          if (score > score_max) then score_max := score;  // v0.4b
        end;
    end;
  // now index must be > 0 in any case (except for H2, or all tags have been cleared)
  if (score_max < 33) then   // v0.4b   fallback for unbranched fragments (use terminal atom!)
    begin
      for i := 1 to ndl_n_atoms do
        begin
          n_nb := ndl_atom^[i].neighbor_count;
          if (n_nb = 1) and (ndl_atom^[i].heavy) and (ndl_atom^[i].tag) then index := i;   
        end;                     // v0.4c (added heavy atom check); v0.4d (added tag check)
    end;
  if index = 0 then inc(index);  // just to be sure...
  find_ndl_ref_atom := index;  
end;


procedure cv_init;  // new in v0.3j
var
  i : integer;
begin
  if (cv = nil) then exit;
  fillchar(cv^,sizeof(connval_type),0);
  for i := 1 to ndl_n_atoms do
    begin
      cv^[i].def := ndl_atom^[i].neighbor_count;
    end;
end;


function cv_count:integer;  // new in v0.3j, modified in v0.3m
var
  i, j : integer;
  cvlist : array[1..max_atoms] of integer;
  cvdef  : integer;
  isnew  : boolean;
  entries : integer;
begin
  if (cv = nil) then 
    begin
      cv_count := 0;
      exit;
    end;
  fillchar(cvlist,sizeof(cvlist),0);
  entries := 0;
  for i := 1 to ndl_n_atoms do
    begin
      if (ndl_atom^[i].heavy = true) then 
        begin
          cvdef := cv^[i].def;
          isnew := true;
          if (entries > 0) then
            begin
              for j := 1 to entries do
                if (cvlist[j] = cvdef) then isnew := false;
            end;
          if isnew then
            begin
              inc(entries);
              cvlist[entries] := cvdef;
            end;
          // now we have a list of unique connection values
        end;
    end;
  cv_count := entries;
end;


function cv_iterate(n_cv_prev: integer): integer;  // new in v0.3j, modified in v0.3m
var
  i, j : integer;
  nb : neighbor_rec;
  nnb : integer;
  nsum : integer;
  n_cv : integer;
begin
  if (cv = nil) or (ndl_n_atoms = 0) then exit;
  // update the connection values (Morgan algorithm)
  for i := 1 to ndl_n_atoms do
    begin
      if (ndl_atom^[i].heavy = true) then 
        begin
          nb   := get_ndl_neighbors(i);
          nnb  := ndl_atom^[i].neighbor_count;
          nsum := 0;
          if (nnb > 0) then
            begin
              for j := 1 to nnb do 
                begin
                  if (ndl_atom^[(nb[j])].heavy = true) then nsum := nsum + cv^[(nb[j])].def;
                end;  
            end;
          cv^[i].tmp := nsum;
        end;
    end;
  n_cv := cv_count;
  if (n_cv > n_cv_prev) then 
    begin
      for i := 1 to ndl_n_atoms do
        begin
          cv^[i].def := cv^[i].tmp;
        end;
    end;    
  cv_iterate := n_cv;
end;


function find_ndl_ref_atom_cv: integer;   // new in v0.3j, modified in v0.3m
var
  i, res, it : integer;
  n_cv, n_cv_prev : integer;
  finished   : boolean;
  cvmax : integer;
begin
  if (ndl_n_atoms = 0) then
    begin
      find_ndl_ref_atom_cv := 0;
      exit;
    end;
  res := 1;
  try
    getmem(cv,sizeof(connval_type));
  except
    on e:Eoutofmemory do
      begin
        res := find_ndl_ref_atom;
        {$IFDEF debug}
        debugoutput('memory allocation for connection values failed, reverting to standard procedure');
        {$ENDIF}
      end;
  end;
  cv_init;
  n_cv_prev := 0;
  it := 0;
  finished := false;
  repeat
    inc(it);  // iteration counter (a safeguard against infinite loops)
    n_cv := cv_iterate(n_cv_prev);
    if (n_cv <= n_cv_prev) then finished := true;
    n_cv_prev := n_cv;
  until (finished or (it > 10000));
  // now that we have canonical connection values (Morgan algorithm),
  // pick the atom with the highest value
  // added in v0.3o: atom must be "tagged"
  cvmax := 0;
  for i := 1 to ndl_n_atoms do
    begin
      //writeln('cv for atom ',i,': ',cv^[i].def);
      if (cv^[i].def > cvmax) and ((ndl_alkene_C(i) = false) or (ez_search = false))
          and (ndl_atom^[i].tag) then   // v0.3o
        begin
          cvmax := cv^[i].def;
          res   := i;
        end;
    end;
  find_ndl_ref_atom_cv := res;  
  try
    if (cv <> nil) then 
      begin
        freemem(cv,sizeof(connval_type));
        cv := nil;
      end;  
  except
    on e:Einvalidpointer do begin end;
  end;
end;


function atomtypes_OK_strict(ndl_a, hst_a : integer):boolean;  // new in v0.2f
var
  ndl_el : str2;
  ndl_atype : str3;
  hst_el : str2;
  hst_atype : str3;
  ndl_nbc : integer;
  hst_nbc : integer;
  ndl_Hexp : smallint;
  ndl_Htot : smallint;  // v0.5a
  hst_Htot : smallint;
  hst_Himp : smallint;  // v0.5a
  res : boolean;
begin
  res := false;   
  ndl_el := ndl_atom^[ndl_a].element;
  ndl_atype := ndl_atom^[ndl_a].atype;
  ndl_nbc := ndl_atom^[ndl_a].neighbor_count;
  ndl_Hexp := ndl_atom^[ndl_a].Hexp;
  ndl_Htot := ndl_atom^[ndl_a].Htot;  // v0.5a
  hst_el := atom^[hst_a].element;
  hst_atype := atom^[hst_a].atype;
  hst_nbc := atom^[hst_a].neighbor_count;
  hst_Htot := atom^[hst_a].Htot;
  hst_Himp := atom^[hst_a].Htot - atom^[hst_a].Hexp;  // v0.5a
  // v0.3o: formal charges must be the same
  if (ndl_atom^[ndl_a].formal_charge <> atom^[hst_a].formal_charge) then 
    begin
      atomtypes_OK_strict := false;
      exit;
    end;
  // v0.3p: isotope nucleon numbers must be the same
  if (ndl_atom^[ndl_a].nucleon_number <> atom^[hst_a].nucleon_number) then 
    begin
      atomtypes_OK_strict := false;
      exit;
    end;
  // v0.3p: radicals numbers must be the same
  if (ndl_atom^[ndl_a].radical_type <> atom^[hst_a].radical_type) then 
    begin
      atomtypes_OK_strict := false;
      exit;
    end;
  if (ndl_atype = hst_atype) then res := true else
    begin
      if (ndl_el = hst_el) then
        begin
          if ((ndl_atom^[ndl_a].arom) and (atom^[hst_a].arom)) then res := true;
          if (ndl_querymol and (ndl_atom^[ndl_a].q_arom) and 
            (atom^[hst_a].arom)) then res := true;   // v0.3p
        end;
    end;
  if (ndl_el = 'A ') and (atom^[hst_a].heavy) then res := true;
  if ndl_el = 'Q ' then
    begin
      if (atom^[hst_a].heavy) and (hst_el <> 'C ') then res := true;
    end;
  if ndl_el = 'X ' then
    begin
      if (hst_el = 'F ') or (hst_el = 'CL') or (hst_el = 'BR') or 
         (hst_el = 'I ') then res := true;
    end;
  // if needle atom has more substituents than haystack atom ==> no match
  if (ndl_nbc > hst_nbc) then res := false;
  // check for explicit hydrogens
  if (ndl_Hexp > hst_Htot) then res := false;
  {$IFDEF debug}  
  if res then debugoutput('atom types OK ('+inttostr(ndl_a)+'/'+inttostr(hst_a)+')')
      else debugoutput('atom types not OK ('+inttostr(ndl_a)+':'+ndl_atype+'/'+inttostr(hst_a)+':'+hst_atype+')');
  {$ENDIF}
  // new in v0.3m: in "fingerprint mode", also query atom symbols must match
  if opt_fp then
    begin
      if ndl_el <> hst_el then res := false;
    end;
  // new in v0.5a: do not allow differing sums of substituents + H(total); this excludes
  // tautomeric forms of (e.g.) N-heteroaromatics
  if ((ndl_nbc + ndl_Htot) <> (hst_nbc + hst_Htot)) then res := false;
  atomtypes_OK_strict := res;
end;

function atomtypes_OK(ndl_a, hst_a : integer):boolean;
var
  ndl_el : str2;
  hst_el : str2;
  ndl_nbc : integer;
  hst_nbc : integer;
  ndl_Hexp : smallint;
  hst_Htot : smallint;
  res : boolean;
begin
  if (ndl_a < 1) or (ndl_a > ndl_n_atoms) or
     (hst_a < 1) or (hst_a > n_atoms) then
    begin
      atomtypes_OK := false;
      exit;
    end;
  // check for opposite charges;  v0.3l, refined in v0.3o, v0.3p
  // except in strict mode, matching pairs of charged+uncharged atoms 
  // are tolerated (this is a feature, not a bug)
  //if ((ndl_atom^[ndl_a].formal_charge <> 0) and (atom^[hst_a].formal_charge <> 0) 
  //and (ndl_atom^[ndl_a].formal_charge <> atom^[hst_a].formal_charge)) then    
  if opt_chg then
    begin
      if (ndl_atom^[ndl_a].formal_charge <> atom^[hst_a].formal_charge) then
        begin
          atomtypes_OK := false;
          exit;
        end;
    end;
  // v0.3p: isotopes must be the same
  if opt_iso then
    begin
      if (ndl_atom^[ndl_a].nucleon_number <> atom^[hst_a].nucleon_number) then
        begin
          atomtypes_OK := false;
          exit;
        end;
    end;
  // v0.3p: radicals must be the same
  if opt_rad then
    begin
      if (ndl_atom^[ndl_a].radical_type <> atom^[hst_a].radical_type) then
        begin
          atomtypes_OK := false;
          exit;
        end;
    end;
  // in exact mode, check if (disconnected) fragment is already tagged; v0.3o
  if opt_exact and (atom^[hst_a].tag = true) then 
    begin
      {$IFDEF debug}
      debugoutput('fragmnet already tagged at '+inttostr(hst_a));
      {$ENDIF}
      atomtypes_OK := false;
      exit;
    end;
  if opt_strict then  // new in v0.2f
    begin
      atomtypes_OK := atomtypes_OK_strict(ndl_a, hst_a);
      exit;
    end;   
  res := false;   
  ndl_el := ndl_atom^[ndl_a].element;
  ndl_nbc := ndl_atom^[ndl_a].neighbor_count;
  ndl_Hexp := ndl_atom^[ndl_a].Hexp;
  hst_el := atom^[hst_a].element;
  hst_nbc := atom^[hst_a].neighbor_count;
  hst_Htot := atom^[hst_a].Htot;
  if (ndl_el = hst_el) then res := true;  // very simplified...
  if (ndl_el = 'A ') and (atom^[hst_a].heavy) then res := true;
  if ndl_el = 'Q ' then
    begin
      if (atom^[hst_a].heavy) and (hst_el <> 'C ') then res := true;
    end;
  if ndl_el = 'X ' then
    begin
      if (hst_el = 'F ') or (hst_el = 'CL') or (hst_el = 'BR') or 
         (hst_el = 'I ') then res := true;
    end;
  // v0.3o: in exact mode, check for identical neighbor_count
  if opt_exact then
    begin
      if (ndl_nbc <> hst_nbc) then 
        begin
          res := false;
          {$IFDEF debug}
          debugoutput('exact match failed: different number of neighbor atoms');
          {$ENDIF}
        end;
    end;
  // if needle atom has more substituents than haystack atom ==> no match
  if (ndl_nbc > hst_nbc) then res := false;
  // check for explicit hydrogens
  if (ndl_Hexp > hst_Htot) then res := false;
  {$IFDEF debug}
  if res then debugoutput('atom types OK ('+inttostr(ndl_a)+'/'+inttostr(hst_a)+')')
      else debugoutput('atom types not OK ('+inttostr(ndl_a)+'/'+inttostr(hst_a)+')');
  {$ENDIF}
  atomtypes_OK := res;
end;

function bondtypes_OK_strict(ndl_b, hst_b : integer):boolean;
var
  ndl_arom  : boolean;
  hst_arom  : boolean;
  ndl_btype : char;
  hst_btype : char;
  ndl_rc    : integer;   // new in v0.3d
  hst_rc    : integer;   // new in v0.3d
  ndl_btopo : shortint;  // new in v0.3d
  res       : boolean;
  {$IFDEF debug}
  na, ha    : string;
  tstr : string;
  {$ENDIF}
begin
  res       := false;
  {$IFDEF debug}
  tstr      := '';  // for debugging purposes only
  {$ENDIF}
  ndl_arom  := ndl_bond^[ndl_b].arom;
  ndl_btype := ndl_bond^[ndl_b].btype;
  ndl_rc    := ndl_bond^[ndl_b].ring_count;
  ndl_btopo := ndl_bond^[ndl_b].topo;
  hst_arom  := bond^[hst_b].arom;
  hst_btype := bond^[hst_b].btype;
  hst_rc    := bond^[hst_b].ring_count;
  {$IFDEF debug}
  if ndl_arom then na := '(ar)' else na := '';
  if hst_arom then ha := '(ar)' else ha := '';
  {$ENDIF}
  if (ndl_arom = true) and (hst_arom = true) then
    begin
      res := true;
      if (ndl_btype = 'T') and (hst_btype <> 'T') then res := false;  // v0.4b
      if (hst_btype = 'T') and (ndl_btype <> 'T') and (ndl_btype <> 'a') then res := false;  // v0.4b
    end;
  if (ndl_arom = false) and (hst_arom = false) then
    begin
      if (ndl_btype = hst_btype) then res := true;
      if (ndl_btype = 'l') and ((hst_btype = 'S') or (hst_btype = 'D')) then res := true;
      if (ndl_btype = 's') and (hst_btype = 'S') then res := true;
      if (ndl_btype = 'd') and (hst_btype = 'D') then res := true;
    end;
  // a little exception:
  if (ndl_arom = false) and (hst_arom = true) then
    begin
      if (ndl_btype = 'A') then res := true;
      if (ndl_btype = 's') or (ndl_btype = 'd') then res := true;
      if (ndl_bond^[ndl_b].q_arom = true) then res := true;  // v0.3p
    end;
  if (ndl_btype = 'a') then res := true;
  // new in v0.3d: strict comparison of topology (and even ring_count!)
  if ((ndl_btopo < btopo_always_any) or (ndl_btopo = btopo_exact_rc)) then
    begin
      if (ndl_rc <> hst_rc) then 
        begin
          res := false;                    // this excludes further ring annulations as well as
          {$IFDEF debug}                   // open-chains query structures to be found in rings
          tstr := ' ringcount mismatch ('+inttostr(ndl_rc)+'/'+inttostr(hst_rc)+')';   
          {$ENDIF}
        end;
    end else
    begin
      if (ndl_btopo = btopo_excess_rc) and (hst_rc <= ndl_rc) then 
        begin
          res := false;
          {$IFDEF debug}
          tstr := ' ringcount mismatch ('+inttostr(ndl_rc)+'/'+inttostr(hst_rc)+')';
          {$ENDIF}
        end;
    end;
  {$IFDEF debug}
  if res then debugoutput('bond types OK ('+inttostr(ndl_b)+':'+ndl_btype+na+'/'+inttostr(hst_b)+':'+hst_btype+ha+')') else
    debugoutput('bond types not OK ('+inttostr(ndl_b)+':'+ndl_btype+na+'/'+inttostr(hst_b)+':'+hst_btype+ha+tstr+')');
  {$ENDIF}
  bondtypes_OK_strict := res;
end;

function bondtypes_OK(ndl_b, hst_b : integer):boolean;
var
  ndl_arom : boolean;
  hst_arom : boolean;
  ndl_btype : char;
  hst_btype : char;
  ndl_rc    : integer;   // new in v0.3d
  hst_rc    : integer;   // new in v0.3d
  ndl_btopo : shortint;  // new in v0.3d
  res : boolean;
  {$IFDEF debug}
  na, ha : string;
  tstr : string;
  q_arom : string;
  {$ENDIF}
  a1, a2 : integer;
  a1_el, a2_el : str2;
begin
  if (ndl_b < 1) or (ndl_b > ndl_n_bonds) or
     (hst_b < 1) or (hst_b > n_bonds) then
     begin
       bondtypes_OK := false;
       exit;
     end;
  if opt_strict then  // new in v0.2f
    begin
      bondtypes_OK := bondtypes_OK_strict(ndl_b, hst_b);
      exit;
    end;
  res := false;
  {$IFDEF debug}
  tstr := '';  // for debug purposes only
  q_arom := 'false';
  {$ENDIF}
  ndl_arom := ndl_bond^[ndl_b].arom;
  ndl_btype := ndl_bond^[ndl_b].btype;
  hst_arom := bond^[hst_b].arom;
  hst_btype := bond^[hst_b].btype;
  ndl_rc    := ndl_bond^[ndl_b].ring_count;
  hst_rc    := bond^[hst_b].ring_count;
  ndl_btopo := ndl_bond^[ndl_b].topo;
  {$IFDEF debug}
  if ndl_arom then na := '(ar)' else na := '';
  if hst_arom then ha := '(ar)' else ha := '';
  if ndl_bond^[ndl_b].q_arom then q_arom := 'true';
  {$ENDIF}
  if (ndl_arom = true) and (hst_arom = true) then
    begin
      res := true;
      if (ndl_btype = 'T') and (hst_btype <> 'T') then res := false;  // v0.4b
      if (hst_btype = 'T') and (ndl_btype <> 'T') and (ndl_btype <> 'a') then res := false;  // v0.4b
    end;
  if (ndl_arom = false) and (hst_arom = false) then
    begin
      if (ndl_btype = hst_btype) then res := true;
      if (ndl_btype = 'l') and ((hst_btype = 'S') or (hst_btype = 'D')) then res := true;
      if (ndl_btype = 's') and (hst_btype = 'S') then res := true;
      if (ndl_btype = 'd') and (hst_btype = 'D') then res := true;
    end;
  // a little exception:
  if (ndl_arom = false) and (hst_arom = true) then
    begin
      if (ndl_btype = 'A') then res := true;
      if (ndl_btype = 's') or (ndl_btype = 'd') then res := true;
      if (ndl_btype = 'D') then   // added in 0.2d: do not accept C=O etc. as C-O/arom
        begin
          a1 := ndl_bond^[ndl_b].a1;
          a2 := ndl_bond^[ndl_b].a2;
          a1_el := ndl_atom^[a1].element;
          a2_el := ndl_atom^[a2].element;
          if not ((a1_el = 'O ') or (a2_el = 'O ') or
                  (a1_el = 'S ') or (a2_el = 'S ') or
                  (a1_el = 'SE') or (a2_el = 'SE') or
                  (a1_el = 'TE') or (a2_el = 'TE')) then res := true;
        end;
      if (ndl_bond^[ndl_b].q_arom = true) then res := true;  // v0.3p
    end;
  if (ndl_btype = 'a') then res := true;
  // new in v0.3d: obey topology requirements in query structure
  if ((ndl_btopo <> btopo_any) and (ndl_btopo <> btopo_always_any)) then
    begin
      if (ndl_btopo = btopo_ring)  and (hst_rc = 0) then res := false;
      if (ndl_btopo = btopo_chain) and (hst_rc > 0) then res := false;
      if (ndl_btopo = btopo_excess_rc) and (hst_rc <= ndl_rc) then res := false;
      if (ndl_btopo = btopo_exact_rc)  and (hst_rc <> ndl_rc) then res := false;
      {$IFDEF debug}
      if res = false then tstr := ' bond topology mismatch '+inttostr(ndl_rc)+'/'+inttostr(hst_rc);
      {$ENDIF}
    end;
  {$IFDEF debug}
  if res then debugoutput('bond types OK ('+inttostr(ndl_b)+':'+ndl_btype+na+'/'
      +inttostr(hst_b)+':'+hst_btype+ha+')') else
    debugoutput('bond types not OK ('+inttostr(ndl_b)+':'+ndl_btype+na+'/'
      +inttostr(hst_b)+':'+hst_btype+ha+tstr+'), q_arom:'+q_arom);
  {$ENDIF}
  bondtypes_OK := res;
end;

function matrix_OK(m:matchmatrix;ndl_dim:integer;hst_dim:integer):boolean;
var   // new, recursive version in v0.2i: can handle up to max_neighbors substituents
  mr : boolean;
  lm : matchmatrix;
  i, ii, j, lndl_dim, lhst_dim : integer;
begin
  if (ndl_dim < 1) or (ndl_dim > max_neighbors) or
     (hst_dim < 1) or (hst_dim > max_neighbors) or
     (ndl_dim > hst_dim) then
    begin
      matrix_OK := false;
      exit;
    end;
  mr := false;
  if (ndl_dim = 1) then
    begin
      for i := 1 to hst_dim do if m[1,i] then mr := true;
    end else
    begin
      for i := 1 to hst_dim do
        begin
          if m[1,i] then
            begin  // write remaining fields into a new matchmatrix which is smaller by 1x1
              fillchar(lm,sizeof(lm),false);
              for j := 2 to ndl_dim do
                begin
                  lhst_dim := 0;
                  for ii := 1 to hst_dim do
                    begin
                      if (ii <> i) then
                        begin
                          inc(lhst_dim);
                          lm[(j-1),lhst_dim] := m[j,ii];
                        end;
                    end;
                end;
              lndl_dim := ndl_dim - 1;  
              if matrix_OK(lm,lndl_dim,lhst_dim) then   // recursive call to matrix_OK
                begin
                  matrix_OK := true;
                  exit;  // stop any further work immediately
                end;  
            end;
        end;
    end;
  matrix_OK := mr;
end;

function is_flat(angle_deg:double):boolean;  // new in v0.3j
begin
 if ((abs(angle_deg) > 5) and (abs(angle_deg) < 175)) then
   is_flat := false else is_flat := true;
end;

function chirality_OK(ndl_cp, hst_cp : chirpath_type):boolean;
var
  res : boolean;
  ndl_ct, hst_ct : double;
  ndl_ct_deg, hst_ct_deg : double;
  np1,np2,np3,np4 : p_3d;
  hp1,hp2,hp3,hp4 : p_3d;
  level, i : integer;
  up, down, updown : boolean;
  ta1,ta2,ta3,ta4,ba1,ba2 : integer;
begin
  if (opt_nochirality) then  // v0.4e
  begin
    chirality_OK := true;
    exit(true);
  end;
  res := true;
  // fill temporary atom variables
  ta1 := ndl_cp[1];  // this is the central atom
  ta2 := ndl_cp[2];
  ta3 := ndl_cp[3];
  ta4 := ndl_cp[4];
  // first, get the central atom of the needle
  np2.x := ndl_atom^[ta1].x;
  np2.y := ndl_atom^[ta1].y;
  np2.z := ndl_atom^[ta1].z;
  // next, do the same for all 3 substituent atoms
  np1.x := ndl_atom^[ta2].x;
  np1.y := ndl_atom^[ta2].y;
  np1.z := ndl_atom^[ta2].z;
  np3.x := ndl_atom^[ta3].x;
  np3.y := ndl_atom^[ta3].y;
  np3.z := ndl_atom^[ta3].z;
  np4.x := ndl_atom^[ta4].x;
  np4.y := ndl_atom^[ta4].y;
  np4.z := ndl_atom^[ta4].z;
  // now check all needle bonds if we should care about up/down bonds
  level  := 0;
  updown := false;
  up     := false;
  down   := false;
  if (ndl_n_bonds > 0) then
    begin
      for i := 1 to ndl_n_bonds do
        begin
          if (ndl_bond^[i].stereo = bstereo_up) or (ndl_bond^[i].stereo = bstereo_down) then
            begin
              ba1 := ndl_bond^[i].a1; ba2 := ndl_bond^[i].a2;
              if (ba1 = ta1) and (ndl_bond^[i].stereo = bstereo_up) then
                begin
                  up := true;
                  if ((ba2 = ta2) or (ba2 = ta3) or (ba2 = ta4)) then
                    begin
                      updown := true;
                      if (ba2 = ta2) then np1.z := np1.z + 0.8;
                      if (ba2 = ta3) then np3.z := np3.z + 0.8;
                      if (ba2 = ta4) then np4.z := np4.z + 0.8;
                    end else level := level + 1;
                end;
              if (ba1 = ta1) and (ndl_bond^[i].stereo = bstereo_down) then
                begin
                  down := true;
                  if ((ba2 = ta2) or (ba2 = ta3) or (ba2 = ta4)) then
                    begin
                      updown := true;
                      if (ba2 = ta2) then np1.z := np1.z - 0.8;
                      if (ba2 = ta3) then np3.z := np3.z - 0.8;
                      if (ba2 = ta4) then np4.z := np4.z - 0.8;
                    end else level := level - 1;
                end;
              if (ba2 = ta1) and (ndl_bond^[i].stereo = bstereo_up) then
                begin
                  down := true;
                  if ((ba1 = ta2) or (ba1 = ta3) or (ba1 = ta4)) then
                    begin
                      updown := true;
                      if (ba1 = ta2) then np1.z := np1.z - 0.8;
                      if (ba1 = ta3) then np3.z := np3.z - 0.8;
                      if (ba1 = ta4) then np4.z := np4.z - 0.8;
                    end else level := level - 1;
                end;
              if (ba2 = ta1) and (ndl_bond^[i].stereo = bstereo_down) then
                begin
                  up := true;
                  if ((ba1 = ta2) or (ba1 = ta3) or (ba1 = ta4)) then
                    begin
                      updown := true;
                      if (ba1 = ta2) then np1.z := np1.z + 0.8;
                      if (ba1 = ta3) then np3.z := np3.z + 0.8;
                      if (ba1 = ta4) then np4.z := np4.z + 0.8;
                    end else level := level + 1;
                end;
            end;  
        end;  // for i ...
      if (updown = false) and (level <> 0) then
        begin
          if (level > 0) then np2.z := np2.z + 0.3;
          if (level < 0) then np2.z := np2.z - 0.3;
        end else
        begin
          if up   then np2.z := np2.z + 0.1;
          if down then np2.z := np2.z - 0.1;
        end;
    end;  
  // fill temporary atom variables again
  ta1 := hst_cp[1];
  ta2 := hst_cp[2];
  ta3 := hst_cp[3];
  ta4 := hst_cp[4];
  // then, get the central atom of the haystack
  hp2.x := atom^[ta1].x;
  hp2.y := atom^[ta1].y;
  hp2.z := atom^[ta1].z;
  // next, do the same for all 3 substituent atoms
  hp1.x := atom^[ta2].x;
  hp1.y := atom^[ta2].y;
  hp1.z := atom^[ta2].z;
  hp3.x := atom^[ta3].x;
  hp3.y := atom^[ta3].y;
  hp3.z := atom^[ta3].z;
  hp4.x := atom^[ta4].x;
  hp4.y := atom^[ta4].y;
  hp4.z := atom^[ta4].z;
  // now check all haystack bonds if we should care about up/down bonds
  level  := 0;
  updown := false;
  up     := false;
  down   := false;
  if (n_bonds > 0) then
    begin
      for i := 1 to n_bonds do
        begin
          if (bond^[i].stereo = bstereo_up) or (bond^[i].stereo = bstereo_down) then
            begin
              ba1 := bond^[i].a1; ba2 := bond^[i].a2;
              if (ba1 = ta1) and (bond^[i].stereo = bstereo_up) then
                begin
                  up := true;
                  if ((ba2 = ta2) or (ba2 = ta3) or (ba2 = ta4)) then
                    begin
                      updown := true;
                      if (ba2 = ta2) then hp1.z := hp1.z + 0.8;
                      if (ba2 = ta3) then hp3.z := hp3.z + 0.8;
                      if (ba2 = ta4) then hp4.z := hp4.z + 0.8;
                    end else level := level + 1;
                end;
              if (ba1 = ta1) and (bond^[i].stereo = bstereo_down) then
                begin
                  down := true;
                  if ((ba2 = ta2) or (ba2 = ta3) or (ba2 = ta4)) then
                    begin
                      updown := true;
                      if (ba2 = ta2) then hp1.z := hp1.z - 0.8;
                      if (ba2 = ta3) then hp3.z := hp3.z - 0.8;
                      if (ba2 = ta4) then hp4.z := hp4.z - 0.8;
                    end else level := level - 1;
                end;
              if (ba2 = ta1) and (bond^[i].stereo = bstereo_up) then
                begin
                  down := true;
                  if ((ba1 = ta2) or (ba1 = ta3) or (ba1 = ta4)) then
                    begin
                      updown := true;
                      if (ba1 = ta2) then hp1.z := hp1.z - 0.8;
                      if (ba1 = ta3) then hp3.z := hp3.z - 0.8;
                      if (ba1 = ta4) then hp4.z := hp4.z - 0.8;
                    end else level := level - 1;
                end;
              if (ba2 = ta1) and (bond^[i].stereo = bstereo_down) then
                begin
                  up := true;
                  if ((ba1 = ta2) or (ba1 = ta3) or (ba1 = ta4)) then
                    begin
                      updown := true;
                      if (ba1 = ta2) then hp1.z := hp1.z + 0.8;
                      if (ba1 = ta3) then hp3.z := hp3.z + 0.8;
                      if (ba1 = ta4) then hp4.z := hp4.z + 0.8;
                    end else level := level + 1;
                end;
            end;  
        end;  // for i ...
      if (updown = false) and (level <> 0) then
        begin
          if (level > 0) then hp2.z := hp2.z + 0.3;
          if (level < 0) then hp2.z := hp2.z - 0.3;
        end else
        begin
          if up   then hp2.z := hp2.z + 0.1;
          if down then hp2.z := hp2.z - 0.1;
        end;
    end;  
  // get the pseudo-torsion angles
  ndl_ct := ctorsion(np1,np2,np3,np4);
  hst_ct := ctorsion(hp1,hp2,hp3,hp4);
  ndl_ct_deg := radtodeg(ndl_ct);
  hst_ct_deg := radtodeg(hst_ct);
  // now do a plausibility check and finally check the sense
  // (clockwise or counterclockwise)
  (*
  if (abs(ndl_ct_deg) > 5) and (abs(ndl_ct_deg) < 175) and
     (abs(hst_ct_deg) > 5) and (abs(hst_ct_deg) < 175) and
     (ndl_ct_deg * hst_ct_deg < 0) then res := false;
  *)
  if (not is_flat(ndl_ct_deg)) and
     (not is_flat(hst_ct_deg)) and
     (ndl_ct_deg * hst_ct_deg < 0) then res := false;
  if rs_strict then
    begin
      if (is_flat(ndl_ct_deg) and (not is_flat(hst_ct_deg))) or
         (is_flat(hst_ct_deg) and (not is_flat(ndl_ct_deg))) or
         (ndl_ct_deg * hst_ct_deg < 0) then res := false;
    end;
  chirality_OK := res;
end;

function ndl_maybe_chiral(na:integer):boolean;  // new in v0.3h
var
  res  : boolean;
  el   : str2;
  at   : str3;
  n_nb : integer;
begin
  res  := false;
  el   := ndl_atom^[na].element;
  at   := ndl_atom^[na].atype;
  n_nb := ndl_atom^[na].neighbor_count;
  if (at = 'C3 ') and (n_nb > 2) then res := true;
  if (el = 'N ') then
    begin
      if (at = 'N3+') and (n_nb = 4) then res := true;
    end;
  if (el = 'S ') then
    begin    // sulfoxide
      if (n_nb = 3) and (ndl_hetatom_count(na) = 1) then res := true;
    end;
  if (el = 'P ') or (el = 'AS') then   // "As" added in v0.3j
    begin
      if (n_nb > 3) then res := true;  // are we missing something here?
      if (ndl_hetatom_count(na) >= 2) then res := false;  // v0.3m; ignore phosphates etc.
    end;
  ndl_maybe_chiral := res;
end;


function is_matching(ndl_xmp:matchpath_type;hst_xmp:matchpath_type):boolean;
var
  i, j, k, l, m : integer;
  ndl_n_nb : integer;
  n_nb : integer;
  ndl_a : integer;
  hst_a : integer;
  ndl_b : integer;
  hst_b : integer;
  prev_ndl_a : integer;
  prev_hst_a : integer;
  next_ndl_a : integer;
  next_hst_a : integer;
  ndl_nb : neighbor_rec;
  hst_nb : neighbor_rec;
  mm : matchmatrix;  
  ndl_mp_len : integer;
  hst_mp_len : integer;
  ndl_mp : matchpath_type;
  hst_mp : matchpath_type;
  emptyline : boolean;
  res : boolean;
  ndl_cis, hst_cis : boolean;
  na1,na2,na3,na4 : integer;  // v0.3d
  ha1,ha2,ha3,ha4 : integer;  // atom variables for E/Z check
  prev_ndl_b : integer;       //
  prev_hst_b : integer;       //
  p1,p2,p3,p4 : p_3d;
  //hst_torsion, ndl_torsion : double;
  ncp,hcp : chirpath_type;
  n_hits : integer;
  n_singlehits : integer;
  {$IFDEF debug}
  tmpstr : string;
  {$ENDIF}
begin
  // initialize local matchpath variables
  fillchar(ndl_mp,sizeof(matchpath_type),0);
  fillchar(hst_mp,sizeof(matchpath_type),0);
  // copy content of external variables into local ones
  for i := 1 to max_matchpath_length do
    begin
      ndl_mp[i] := ndl_xmp[i];
      hst_mp[i] := hst_xmp[i];
    end;
  ndl_mp_len := matchpath_length(ndl_mp);
  hst_mp_len := matchpath_length(hst_mp);
  if ndl_mp_len <> hst_mp_len then
    begin
      // this should never happen....
      {$IFDEF debug}
      debugoutput('needle and haystack matchpaths are of different length');
      {$ENDIF}
      is_matching := false;
      exit;
    end;  
  ndl_a := ndl_mp[ndl_mp_len];
  hst_a := hst_mp[hst_mp_len];
  ndl_atom^[ndl_a].tag := false; // new in v0.3o: mark the last needle atom as "visited"
  ndl_b := 0;
  hst_b := 0;
  prev_ndl_a := 0;
  prev_hst_a := 0;
  if ndl_mp_len > 1 then
    begin
      prev_ndl_a := ndl_mp[ndl_mp_len-1];
      prev_hst_a := hst_mp[hst_mp_len-1];
    end;
  // if geometry checking is on, check it here
  if (ez_search = true) and (ndl_mp_len > 3) then
    begin
      na1 := ndl_mp[ndl_mp_len];
      na2 := ndl_mp[ndl_mp_len-1];
      na3 := ndl_mp[ndl_mp_len-2];
      na4 := ndl_mp[ndl_mp_len-3];
      ha1 := hst_mp[hst_mp_len];
      ha2 := hst_mp[hst_mp_len-1];
      ha3 := hst_mp[hst_mp_len-2];
      ha4 := hst_mp[hst_mp_len-3];
      prev_ndl_b := get_ndl_bond(na2,na3);
      prev_hst_b := get_bond(ha2,ha3);
      if (ndl_bond^[prev_ndl_b].btype = 'D') and
         (bond^[prev_hst_b].arom = false) and
         ((atom^[ha2].element = 'C ') or (atom^[ha2].element = 'N ')) and
         ((atom^[ha3].element = 'C ') or (atom^[ha3].element = 'N ')) then  // v0.3g; check C=C, C=N, N=N bonds
        begin
          p1.x := atom^[ha1].x; p1.y := atom^[ha1].y; p1.z := atom^[ha1].z;
          p2.x := atom^[ha2].x; p2.y := atom^[ha2].y; p2.z := atom^[ha2].z;
          p3.x := atom^[ha3].x; p3.y := atom^[ha3].y; p3.z := atom^[ha3].z;
          p4.x := atom^[ha4].x; p4.y := atom^[ha4].y; p4.z := atom^[ha4].z;
          hst_cis := is_cis(p1,p2,p3,p4);
          //hst_torsion := torsion(p1,p2,p3,p4);
          p1.x := ndl_atom^[na1].x; p1.y := ndl_atom^[na1].y; p1.z := ndl_atom^[na1].z;
          p2.x := ndl_atom^[na2].x; p2.y := ndl_atom^[na2].y; p2.z := ndl_atom^[na2].z;
          p3.x := ndl_atom^[na3].x; p3.y := ndl_atom^[na3].y; p3.z := ndl_atom^[na3].z;
          p4.x := ndl_atom^[na4].x; p4.y := ndl_atom^[na4].y; p4.z := ndl_atom^[na4].z;
          //ndl_torsion := torsion(p1,p2,p3,p4);
          ndl_cis := is_cis(p1,p2,p3,p4);
          if ndl_cis <> hst_cis then
            begin
              {$IFDEF debug}
              debugoutput('E/Z geometry mismatch');
              {$ENDIF}
              is_matching := false;
              exit;
            end;          
        end;
    end;  // end of E/Z geometry check
  // check whatever can be checked as early as now:
  // e.g. different elements or more substituents on needle atom than on haystack
  if (not atomtypes_OK(ndl_a, hst_a)) then
    begin
      is_matching := false;
      exit;
    end;
  // positive scenarios, e.g. one-atom fragments  (v0.3o)
  if (atom^[hst_a].neighbor_count = 0) and (ndl_atom^[ndl_a].neighbor_count = 0) then
    begin
      if atomtypes_OK(ndl_a, hst_a) then
        begin
          is_matching := true;
          atom^[hst_a].tag := true;
          if (use_gmm and valid_gmm) then gmm^[ndl_a,hst_a] := true;  // v0.4a
        end else is_matching := false;
        exit;
    end; 
  // and other possibilities:
  ndl_b := get_ndl_bond(prev_ndl_a,ndl_a);
  hst_b := get_bond(prev_hst_a,hst_a);
  {$IFDEF debug}
  debugoutput('Now checking atoms '+inttostr(ndl_a)+'/'+inttostr(hst_a)+', bonds '+inttostr(ndl_b)+'/'+inttostr(hst_b));
  {$ENDIF}
  if (ndl_b > 0) and (hst_b > 0) then
    begin
      // do a quick check if bond types match
      if (not bondtypes_OK(ndl_b, hst_b)) then
        begin
          {$IFDEF debug}
          debugoutput('  failed match of bonds '+inttostr(ndl_b)+'/'+inttostr(hst_b));
          {$ENDIF}
          is_matching := false;
          exit;
        end;
    end;
  // v0.4b: reversed checks a) and a.1), see below
  // a.1) haystack fragment forms a ring, but needle does not;  v0.3m
  if (matchpath_pos(ndl_a,ndl_mp) = matchpath_length(ndl_mp)) and
     (matchpath_pos(hst_a,hst_mp) < matchpath_length(hst_mp)) then
     begin
       is_matching := false;
       {$IFDEF debug}
       debugoutput('  haystack forms a ring and needle does not at '+inttostr(hst_a));
       {$ENDIF}
       exit;
     end; 
  // a) we reached the end of our needle fragment (and atom/bond types match)
  if (ndl_atom^[ndl_a].neighbor_count = 1) and (atomtypes_OK(ndl_a, hst_a)) and
     (bondtypes_OK(ndl_b, hst_b)) then
     begin
       is_matching := true;
       {$IFDEF debug}
       debugoutput('  ==> end of needle fragment at atom '+inttostr(ndl_a)+' (match)');
       {$ENDIF}
       if (use_gmm and valid_gmm) then gmm^[ndl_a,hst_a] := true;  // v0.4a
       exit;
     end;    
  // b) a ring is formed (ndl_a is already in the path) and atom/bond types match
  if (matchpath_pos(ndl_a,ndl_mp) > 0) and (matchpath_pos(ndl_a,ndl_mp) < matchpath_length(ndl_mp)) then
    begin
      if (matchpath_pos(ndl_a,ndl_mp) = matchpath_pos(hst_a, hst_mp)) and
         (atomtypes_OK(ndl_a, hst_a)) and
         (bondtypes_OK(ndl_b, hst_b)) then      
         begin
           // 1st chirality check
           if (matchpath_pos(ndl_a,ndl_mp) > 1) and (rs_search or ndl_atom^[ndl_a].stereo_care)
           and ndl_maybe_chiral(ndl_a) then                // new in v0.3h
             begin
               na1 := ndl_a;  // the (potential) chiral center (v0.3f)
               na2 := ndl_mp[(matchpath_pos(ndl_a,ndl_mp)-1)];
               na3 := ndl_mp[(matchpath_pos(ndl_a,ndl_mp)+1)];
               na4 := ndl_mp[(matchpath_length(ndl_mp)-1)];
               ha1 := hst_a;  
               ha2 := hst_mp[(matchpath_pos(hst_a,hst_mp)-1)];
               ha3 := hst_mp[(matchpath_pos(hst_a,hst_mp)+1)];
               ha4 := hst_mp[(matchpath_length(hst_mp)-1)];
               fillchar(ncp,sizeof(chirpath_type),0);
               fillchar(hcp,sizeof(chirpath_type),0);
               ncp[1] := na1; ncp[2] := na2; ncp[3] := na3; ncp[4] := na4;
               hcp[1] := ha1; hcp[2] := ha2; hcp[3] := ha3; hcp[4] := ha4;
               if not chirality_OK(ncp,hcp) then
                 begin
                   {$IFDEF debug}
                   debugoutput('chirality check failed at ring junction');
                   {$ENDIF}
                   is_matching := false;
                   exit;
                 end else
                 begin
                   {$IFDEF debug}
                   debugoutput('chirality check succeeded at ring junction');
                   {$ENDIF}
                 end;
             end;  // end of 1st chirality check
           is_matching := true;
           if (use_gmm and valid_gmm) then gmm^[ndl_a,hst_a] := true;  // v0.4a
           {$IFDEF debug}
           debugoutput('matchpath forms ring at: '+inttostr(ndl_a)+' (match)');
           {$ENDIF}
           exit;
         end else
         begin
           is_matching := false;
           {$IFDEF debug}
           debugoutput('matchpath forms ring at: '+inttostr(ndl_a)+' (no match)');
           {$ENDIF}
           exit;
         end;
    end;
  // in all other cases, do the hard work:
  // first, get all heavy-atom neighbors of needle and haystack;
  // at the beginning of the search, this means all neighbors, then it means
  // all but the previous atom (where we came from)
  fillchar(ndl_nb, sizeof(neighbor_rec),0);
  fillchar(hst_nb, sizeof(neighbor_rec),0);
  if (matchpath_length(ndl_mp) = 1) then
    begin
      ndl_n_nb := ndl_atom^[ndl_a].neighbor_count;
      n_nb := atom^[hst_a].neighbor_count;
      ndl_nb := get_ndl_neighbors(ndl_a);
      hst_nb := get_neighbors(hst_a);
    end else
    begin
      ndl_n_nb := ndl_atom^[ndl_a].neighbor_count -1;
      n_nb := atom^[hst_a].neighbor_count -1;
      ndl_nb := get_ndl_nextneighbors(ndl_a,prev_ndl_a);
      hst_nb := get_nextneighbors(hst_a,prev_hst_a);
    end;
  // v0.3o: mark all neighbor atoms as "visited"
  for i := 1 to ndl_n_nb do
    ndl_atom^[(ndl_nb[i])].tag := false;
  // now that the neighbor-arrays are filled, get all
  // combinations of matches recursively;
  // first, initialize the match matrix
  fillchar(mm,sizeof(mm),false);    // new in v0.2i
  // make sure there are not too many neighbors (max. max_neighbors)  
  if (ndl_n_nb > max_neighbors) or (n_nb > max_neighbors) then       // updated in v0.2i
    begin
      {$IFDEF debug}
      debugoutput('too many neighbors - exiting');
      {$ENDIF}
      is_matching := false;
      exit;
    end;
  // check if matchpath is not already filled up
  if matchpath_length(ndl_mp) = max_matchpath_length then
    begin
      {$IFDEF debug}
      debugoutput('matchpath too long - exiting');
      {$ENDIF}
      is_matching := false;
      exit;
    end;
  // next, check which chain of the needle matches which chain of the haystack 
  for i := 1 to ndl_n_nb do
    begin
      emptyline := true;
      next_ndl_a := ndl_nb[i];
      for j := 1 to n_nb do
        begin
          next_hst_a := hst_nb[j];
          ndl_mp[ndl_mp_len+1] := next_ndl_a;
          hst_mp[hst_mp_len+1] := next_hst_a;
          if is_matching(ndl_mp,hst_mp) then  // recursive function call
            begin
              inc(recursion_depth);     // new in v0.3p: limit the recursion depth
              if (max_match_recursion_depth > 0) and 
              (recursion_depth > max_match_recursion_depth) then
                begin
                  if (opt_verbose) then
                    writeln('Warning: max. number of match recursions (',
                      max_match_recursion_depth,
                      ') reached, reverting to non-exhaustive match');
                  is_matching := true;
                  if (use_gmm and valid_gmm) then gmm^[ndl_a,hst_a] := true;  // v0.4a
                  exit;
                end;                    // end of v0.3p recursion depth checking 
              mm[i,j] := true; 
              emptyline := false;
            end;
        end;
      // if a needle substituent does not match any of the haystack substituents,
      // stop any further work immediately
      if emptyline then 
        begin
          is_matching := false;
          exit;
        end;  
    end; 
  // finally, check the content of the matrix
  res := matrix_OK(mm,ndl_n_nb,n_nb);
  // optional: chirality check
  if (res and (rs_search or (use_gmm and valid_gmm) or ndl_atom^[ndl_a].stereo_care) 
    and ndl_maybe_chiral(ndl_a)) then // v0.4a
    begin
      // first, we have to clean up the match matrix in order to remove
      // "impossible" multiple matches (new in v0.3h)
      for i := 1 to 3 do
        begin
          for j := 1 to max_neighbors do  // haystack dimension
            begin
              n_hits := 0;
              l    := 0;
              for k := 1 to max_neighbors do  // needle dimension
                begin
                  if mm[k,j] then
                    begin
                      inc(n_hits);
                      l := k;
                    end;
                end;
              if (n_hits = 1) then  // a unique match ==> kick out any other match at this pos.
                begin
                  for m := 1 to max_neighbors do
                    begin
                      if (m <> j) then
                        begin
                          if mm[l,m] then   // v0.4a
                            begin
                              if (use_gmm and valid_gmm) then gmm^[ndl_nb[l],hst_nb[m]] := false;              
                            end;
                          mm[l,m] := false;
                        end;
                    end;
                end;
            end;
        end;
      // end of match matrix clean-up
      if (prev_ndl_a > 0) then
        begin
          n_singlehits := 1;
          ncp[2] := prev_ndl_a;
          hcp[2] := prev_hst_a;
        end else
        begin
          n_singlehits := 0;
        end;
      ncp[1] := ndl_a;
      hcp[1] := hst_a;  
      i := 0; l := 0;
      while (n_singlehits < 3) and (i < 4) do
        begin
          inc(i);
          n_hits := 0;
          for k := 1 to n_nb do
            begin
              if mm[i,k] then 
                begin
                  inc(n_hits);
                  l := k;
                end;
            end;
          if (n_hits = 1) then 
            begin
              inc(n_singlehits);
              ncp[(n_singlehits+1)] := ndl_nb[i];  
              hcp[(n_singlehits+1)] := hst_nb[l];
            end;  
        end;
      if (n_singlehits = 3) then
        begin
          if not chirality_OK(ncp,hcp) then
            begin
              {$IFDEF debug}
              debugoutput('chirality check failed');
              {$ENDIF}
              res := false;
            end else
            begin
              {$IFDEF debug}
              debugoutput('chirality check OK');
              {$ENDIF}
            end;
        end;
    end;
  is_matching := res;
  {$IFDEF debug}
  if res then tmpstr := ' MATCH' else tmpstr := ' NO MATCH';
  debugoutput('result for atoms '+inttostr(ndl_a)+'/'+inttostr(hst_a)+', bonds '+inttostr(ndl_b)+'/'+inttostr(hst_b)+':'+tmpstr);
  {$ENDIF}
  if (res and use_gmm and valid_gmm) then gmm^[ndl_a,hst_a] := true;  // v0.4a
end;


function quick_match:boolean;  // added in v0.2c
var
  res : boolean;
  i : integer;
  ndl_atype : str3;
  ndl_el    : str2;    // v0.3l
  ndl_chg   : integer; // v0.3l
  ndl_rad   : integer; // v0.3p
  ndl_iso   : integer; // v0.3p
begin
  if (ez_search or rs_search) and (ndl_n_heavyatoms > 3) then  // v0.3f, v0.3m, v0.3o
    begin
      quick_match := false;
      exit;
    end;
  if ((ndl_n_atoms < 1) or (n_atoms < 1)) or
     ((ndl_n_atoms > n_atoms) or (ndl_n_bonds > n_bonds)) then
    begin   // just to be sure...
      quick_match := false;
      {$IFDEF debug}
      debugoutput(' ==> quick_match failed');
      {$ENDIF}
      exit;
    end;
  res := true;
  for i := 1 to ndl_n_atoms do
    begin
      //if atom^[i].atype <> ndl_atom^[i].atype then res := false;    // changed in
      if atom^[i].element <> ndl_atom^[i].element then res := false;  // v0.2k
      if opt_chg and (atom^[i].formal_charge <> ndl_atom^[i].formal_charge) then res := false;  // v0.3o, v0.3p
      if opt_iso and (atom^[i].nucleon_number <> ndl_atom^[i].nucleon_number) then res := false;  // v0.3p
      if opt_rad and (atom^[i].radical_type <> ndl_atom^[i].radical_type) then res := false;  // v0.3p
    end;
  {$IFDEF debug}
  if res then debugoutput(' ==> quick_match: atoms OK') else debugoutput(' ==> quick_match: atoms not OK');
  {$ENDIF}
  if ndl_n_bonds > 0 then
    begin
      for i := 1 to ndl_n_bonds do
        begin
          if (ndl_bond^[i].a1 <> bond^[i].a1) or
             (ndl_bond^[i].a2 <> bond^[i].a2) or
             (ndl_bond^[i].btype <> bond^[i].btype) then res := false;
        end;
    end;
  {$IFDEF debug}
  if res then debugoutput(' ==> quick_match: bonds OK') else debugoutput(' ==> quick_match: bonds not OK');
  {$ENDIF}
  // added in v0.2d: special case: needle contains only one heavy atom; refined in v0.3l, v0.3o
  if (ndl_n_heavyatoms = 1) then
    begin
      res := false;  // v0.3p
      // first, find out the element and atom type of the only heavy atom      
      for i := 1 to ndl_n_atoms do
        if ndl_atom^[i].heavy then 
          begin
            ndl_atype := ndl_atom^[i].atype;
            ndl_el    := ndl_atom^[i].element;  // v0.3l
            ndl_chg   := ndl_atom^[i].formal_charge;  // v0.3l
            ndl_iso   := ndl_atom^[i].nucleon_number; // v0.3p
            ndl_rad   := ndl_atom^[i].radical_type;   // v0.3p
          end;
      for i := 1 to n_atoms do
        begin
          if (atom^[i].atype = ndl_atype) and (atom^[i].element = ndl_el) and
             (ndl_chg = atom^[i].formal_charge) then res := true;  // v0.3l, v0.3o
          if res and opt_iso then if atom^[i].nucleon_number <> ndl_iso then res := false; // v0.3p
          if res and opt_rad then if atom^[i].radical_type <> ndl_rad then res := false; // v0.3p
        end;
    end;
  {$IFDEF debug}
  if res then debugoutput(' ==> quick_match succeeded') else debugoutput(' ==> quick_match failed (2)');
  {$ENDIF}
  if (res and use_gmm and valid_gmm and opt_matchnum1) then     // v0.4a
    begin
      for i := 1 to ndl_n_atoms do gmm^[i,i] := true;
    end;
  quick_match := res;
end;

function gmm_collision:boolean; forward;  // v0.4b

procedure perform_match;
var
  i, j : integer;
  //ndl_ref_atom : integer;  // since v0.3j as a global variable
  ndl_n_nb : integer;
  ndl_n_hc : integer;
  n_nb : integer;
  n_hc : integer;
  qm   : boolean;  // v0.3l
begin
  // check for NoStruct (0 atoms);  v0.3l
  if (n_atoms = 0) or (ndl_n_atoms = 0) then
    begin
      matchresult := false;
      {$IFDEF debug}
      debugoutput('NoStruct encountered - aborted match routine');
      {$ENDIF}
      exit;
    end;
  // if we perform an exact match, needle and haystack must have
  // the same number of atoms, bonds, and rings
  if (opt_exact and opt_iso) then
    begin
      if (n_heavyatoms <> ndl_n_heavyatoms) or (n_heavybonds <> ndl_n_heavybonds) then 
        begin
          matchresult := false;
          {$IFDEF debug}
          debugoutput('different number of heavy atoms and/or bonds');
          {$ENDIF}
          exit;
        end;
    end;
  if (opt_exact and (opt_iso = false)) and (n_trueheavyatoms <> ndl_n_trueheavyatoms) then   // v0.4b
    begin
      matchresult := false;
      {$IFDEF debug}
      debugoutput('different number of (true) heavy atoms');
      {$ENDIF}
      exit;
    end;    
  // have a quick look if needle and haystack are identical molfiles
  if ((use_gmm = false) or (valid_gmm = false) or (opt_matchnum = false)) then    // v0.4a
    begin  
      qm := quick_match;  // v0.3l
      if qm then
        begin
          matchresult := true;
          clear_ndl_atom_tags;  // v0.3o
          exit;
        end;
      end;
  // if we have only one heavy atom and quick_match fails, return "false";  v0.3l
  if (qm = false) and (ndl_n_heavyatoms = 1) then
    begin
      matchresult := false;
      exit;
    end;
  {$IFDEF debug}
  debugoutput('needle reference atom: '+inttostr(ndl_ref_atom)+' ('+ndl_atom^[ndl_ref_atom].atype+')');
  {$ENDIF}
  ndl_n_nb := ndl_atom^[ndl_ref_atom].neighbor_count;
  ndl_n_hc := ndl_hetatom_count(ndl_ref_atom);
  {$IFDEF debug}
  debugoutput('neighbor atoms: '+inttostr(ndl_n_nb)+'  heteroatom neighbors: '+inttostr(ndl_n_hc));
  {$ENDIF}
  matchresult := false;
  for j := 1 to max_matchpath_length do
    begin
      ndl_matchpath[j] := 0;
      hst_matchpath[j] := 0;
    end;
  ndl_matchpath[1] := ndl_ref_atom;
  i := 0;
  found_untagged := false;
  while (i < n_atoms) and (matchresult = false) do
    begin
      inc(i);
      recursion_depth := 0;  // v0.3p
      n_nb := atom^[i].neighbor_count;
      n_hc := hetatom_count(i);
      if ((n_nb >= ndl_n_nb) and (n_hc >= ndl_n_hc))
         and (not (use_gmm and valid_gmm and tmp_tag[i])) then     // v0.4a
        begin
          if (use_gmm and valid_gmm) then found_untagged := true;   // v0.4a
          {$IFDEF debug}
          debugoutput('trying atom '+inttostr(i)+'; neighbor atoms: '+inttostr(n_nb)+' heteroatom neighbors: '+inttostr(n_hc));
          {$ENDIF}
          hst_matchpath[1] := i;
          matchresult := is_matching(ndl_matchpath, hst_matchpath);
          if (use_gmm and valid_gmm) then   // v0.4b
            begin
             if gmm_collision then
               begin
                 matchresult := false;
                 clear_gmm;
               end;
            end;
          {$IFDEF debug}
          if matchresult then debugoutput('matching atom in haystack: '+inttostr(i)+' ('+atom^[i].atype+')');
          {$ENDIF}
          if matchresult then atom^[i].tag := true;  // v0.3o; mark this fragment as matched
          tmp_tag[i] := true;    // v0.4a
        end;
    end;
end;


procedure clear_rings;
var
  i : integer;
begin
  n_rings := 0;
  fillchar(ring^, sizeof(ringlist),0);
  for i := 1 to max_rings do  // new in v0.3
    begin
      ringprop^[i].size     := 0;
      ringprop^[i].arom     := false;
      ringprop^[i].envelope := false;
    end;
  if n_atoms > 0 then
    begin
      for i := 1 to n_atoms do atom^[i].ring_count := 0;
    end;
  if n_bonds > 0 then
    begin
      for i := 1 to n_bonds do bond^[i].ring_count := 0;
    end;
end;


function ring_lastpos(s:ringpath_type):integer;
var
  i, rc, rlp : integer;
begin
  rlp := 0;
  if n_rings > 0 then
    begin
      for i := 1 to n_rings do
        begin
          rc := ringcompare(s, ring^[i]);
          if rc_identical(rc) then rlp := i;
        end;
    end;
  ring_lastpos := rlp;
end;


procedure remove_redundant_rings;
var
  i, j, k, rlp : integer;
  tmp_path : ringpath_type;
begin
  if n_rings < 2 then exit;
  for i := 1 to (n_rings - 1) do
    begin
      tmp_path := ring^[i];
      rlp := ring_lastpos(tmp_path);
      while rlp > i do
        begin
          {$IFDEF debug}
          debugoutput('removing redundant ring: '+inttostr(rlp)+' (identical to ring '+inttostr(i)+')');
          {$ENDIF}
          for j := rlp to (n_rings - 1) do
            begin
              ring^[j] := ring^[(j+1)];
              ringprop^[j].size := ringprop^[(j+1)].size;  // new in v0.3
              ringprop^[j].arom := ringprop^[(j+1)].arom;
              ringprop^[j].envelope := ringprop^[(j+1)].envelope;
            end;
          for k := 1 to max_ringsize do ring^[n_rings,k] := 0;
          dec(n_rings);
          rlp := ring_lastpos(tmp_path);
        end;
    end;
end;

function count_aromatic_rings:integer;
var
  i, n : integer;
begin
  n := 0;
    if n_rings > 0 then
      begin
        for i := 1 to n_rings do
          if ringprop^[i].arom then inc(n);
      end;
  count_aromatic_rings := n;
end;


procedure chk_envelopes;  // new in v0.3d
// checks if a ring completely contains one or more other rings
var
  a,i,j,k,l,pl,pli : integer;
  found_atom, found_all_atoms, found_ring : boolean;
begin
  if n_rings < 2 then exit;
  for i := 2 to n_rings do
    begin
      found_ring := false;
      j := 0;
      pli := ringprop^[i].size;  // path_length(ring^[i]);
      while ((j < (i-1)) and (found_ring = false)) do
        begin
          inc(j);
          found_all_atoms := true;
          pl := ringprop^[j].size;  // path_length(ring^[j]);
          for k := 1 to pl do
            begin
              found_atom := false;
              a := ring^[j,k];
              for l := 1 to pli do
                begin
                  if ring^[i,l] = a then found_atom := true;
                end;
              if found_atom = false then found_all_atoms := false;
            end;
          if found_all_atoms then found_ring := true;
        end;
      if found_ring then ringprop^[i].envelope := true;
    end;
end;

procedure update_ringcount;
var
  i, j, a1, a2, b, pl : integer;
begin
  if n_rings > 0 then
    begin
      chk_envelopes;
      for i := 1 to n_rings do
        begin
          if (ringprop^[i].envelope = false) then
            begin
              pl := ringprop^[i].size;  // path_length(ring^[i]);  // v0.3d
              a2 := ring^[i,pl];
              for j := 1 to pl do
                begin
                  a1 := ring^[i,j];
                  inc(atom^[a1].ring_count);
                  b := get_bond(a1,a2);
                  inc(bond^[b].ring_count);
                  a2 := a1;
                end;
            end;
        end;
    end;
end;


function normalize_ionic_bonds:boolean;  // v0.3k
var                                      // changed from a procedure into a function in v0.3m
  i : integer;
  a1, a2 : integer;
  fc1, fc2 : integer;
  bt : char;
  res : boolean; // v0.3m
begin
  res := false; // v0.3m
  if (n_bonds = 0) then 
    begin
      normalize_ionic_bonds := res;
      exit;
    end;
  for i := 1 to n_bonds do
    begin
      a1 := bond^[i].a1;
      a2 := bond^[i].a2;
      bt := bond^[i].btype;
      fc1 := atom^[a1].formal_charge;
      fc2 := atom^[a2].formal_charge;
      if (fc1 * fc2 = -1) and ((bt = 'S') or (bt = 'D')) then
        begin
          atom^[a1].formal_charge := 0;
          atom^[a2].formal_charge := 0;
          if (atom^[a1].atype = 'N3+') then atom^[a1].atype := 'N3 ';  // v0.3m
          if (atom^[a2].atype = 'N3+') then atom^[a2].atype := 'N3 ';  // v0.3m
          if (bt = 'D') then 
            begin   // v0.5a restrict changes from D to T to only C=C, N=N and C=N bonds
              if ((atom^[a1].element = 'C ') and ((atom^[a2].element = 'C ') or (atom^[a2].element = 'N '))) or
                 ((atom^[a1].element = 'N ') and ((atom^[a2].element = 'C ') or (atom^[a2].element = 'N ')))
              then bond^[i].btype := 'T';
            end;
          if (bt = 'S') then bond^[i].btype := 'D';
          res := true;  // v0.3m
        end;
    end;
  normalize_ionic_bonds := res;  // v0.3m (return true if any change was made
end;

procedure chk_wildcard_rings;  // new in v0.3p
// checks if there are any wildcard atom types or bond types
// in a ring of the needle; if yes ==> set the q_arom flag in the
// atom and bond record of all ring members in order to perform the 
// match a bit more generously
var
  i, j, rs : integer;
  a1, a2, b : integer;
  wcr : boolean;
  at : str3;
  bt : char;
begin
  if (ndl_querymol = false) then exit;
  if (ndl_n_rings = 0) then exit;
  // now look for any not-yet-aromatic rings which contain a wildcard
  for i := 1 to ndl_n_rings do
    begin
      wcr := false;
      if (ndl_ringprop^[i].arom = false) then
        begin
          rs := ndl_ringprop^[i].size;
          a2 := ndl_ring^[i,rs];
          for j := 1 to rs do
            begin
              a1 := ndl_ring^[i,j];
              b := get_ndl_bond(a1,a2);
              at := ndl_atom^[a1].atype;
              bt := ndl_bond^[b].btype;
              if (at = 'A  ') or (at = 'Q  ') then wcr := true;
              if (bt = 'l') or (bt = 's') or (bt = 'd') or (bt = 'a') then wcr := true;
              a2 := a1;   
            end;
          if wcr then 
            begin  // if yes, flag all atoms and bonds in this ring as "potentially" aromatic
              {$IFDEF debug}
              debugoutput('wildcard ring found');
              {$ENDIF}
              a2 := ndl_ring^[i,rs];
              for j := 1 to rs do
                begin
                  a1 := ndl_ring^[i,j];
                  b := get_ndl_bond(a1,a2);
                  at := ndl_atom^[a1].atype;
                  bt := ndl_bond^[b].btype;
                  ndl_atom^[a1].q_arom := true;
                  ndl_bond^[b].q_arom := true;
                  a2 := a1;   
                end;
            end;
        end;
    end;
  // and now undo this flagging for all rings which contain no wildcard
  for i := 1 to ndl_n_rings do
    begin
      wcr := false;
      rs := ndl_ringprop^[i].size;
      a2 := ndl_ring^[i,rs];
      for j := 1 to rs do
        begin
          a1 := ndl_ring^[i,j];
          b := get_ndl_bond(a1,a2);
          at := ndl_atom^[a1].atype;
          bt := ndl_bond^[b].btype;
          if (at = 'A  ') or (at = 'Q  ') then wcr := true;
          if (bt = 'l') or (bt = 's') or (bt = 'd') or (bt = 'a') then wcr := true;
          a2 := a1;   
        end;
      if (wcr = false) then 
        begin  // if yes, unflag all atoms and bonds in this ring
          a2 := ndl_ring^[i,rs];
          for j := 1 to rs do
            begin
              a1 := ndl_ring^[i,j];
              b := get_ndl_bond(a1,a2);
              at := ndl_atom^[a1].atype;
              bt := ndl_bond^[b].btype;
              ndl_atom^[a1].q_arom := false;
              ndl_bond^[b].q_arom := false;
              a2 := a1;   
            end;
        end;
    end;
  // some further refinement would be necessary here in order to unflag everything
  // which contains a wildcard but which definitely cannot be aromatic
end;

procedure revert_frag(var frag:frag_rec);
var
  i : integer;
  a_tmp : integer;
  b_tmp : char;
begin
  if frag.size < 2 then exit;
  if (frag.ring = false) then
    begin
      for i := 1 to (frag.size div 2) do
        begin
          a_tmp := frag.a_atomicnum[i];
          b_tmp :=  frag.b_code[i];
          frag.a_atomicnum[i] := frag.a_atomicnum[frag.size-(i-1)];
          frag.a_atomicnum[frag.size-(i-1)] := a_tmp;
          frag.b_code[i] := frag.b_code[frag.size-i];
          frag.b_code[frag.size-i] := b_tmp;
        end;
    end;
end;

procedure order_frag(var frag:frag_rec);
var
  i : integer;
  ref : integer;
begin
  if frag.size < 2 then exit;
  if (frag.ring = false) then
    begin
      ref := 0;
      for i := 1 to (frag.size div 2) do
        begin
          ref := ref + frag.a_atomicnum[i];
          ref := ref - frag.a_atomicnum[frag.size-(i-1)];
          if (ref <> 0) then break;
        end;
      if (ref = 0) then
        begin
          for i := 1 to (frag.size div 2) do
            begin
              ref := ref + ord(frag.b_code[i]);
              ref := ref - ord(frag.b_code[frag.size-i]);
              if (ref <> 0) then break;
            end;
        end;
      if (ref < 0) then revert_frag(frag);
    end;
end;


function mk_fragstr(frag:frag_rec):fragstr;
var
  i : integer;
  res, tmpstr : string;
begin
  res := '';
  tmpstr := '';
  for i := 1 to frag.size do
    begin
      str(frag.a_atomicnum[i],tmpstr);
      res := res + tmpstr;
      tmpstr := '';
      if (i < frag.size) or (frag.ring = true) then tmpstr := frag.b_code[i];
      res := res + tmpstr;
    end;
  mk_fragstr := res;  
end;

// various hash functions have been taken from here:
(**************************************************************************)
(*                                                                        *)
(*          General Purpose Hash Function Algorithms Library              *)
(*                                                                        *)
(* Author: Arash Partow - 2002                                            *)
(* URL: http://www.partow.net                                             *)
(* URL: http://www.partow.net/programming/hashfunctions/index.html        *)
(*                                                                        *)
(* Copyright notice:                                                      *)
(* Free use of the General Purpose Hash Function Algorithms Library is    *)
(* permitted under the guidelines and in accordance with the most current *)
(* version of the Common Public License.                                  *)
(* http://www.opensource.org/licenses/cpl.php                             *)
(*                                                                        *)
(**************************************************************************)


function FNVHash(const str : string) : cardinal;
const 
  FNVPrime = $811C9DC5;
var
  i, res : cardinal;
begin
  res := 0;
  for i := 1 to length(str) do
    begin
      res := res * FNVPrime;
      res := res xor ord(str[i]);
    end;
  FNVHash := res;
end;


function BKDRHash(const str : string) : cardinal;
const 
  seed = 131; (* 31 131 1313 13131 131313 etc... *)
var
  i, res : cardinal;
begin
  res := 0;
  for i := 1 to length(str) do
    begin
      res := (res * seed) + ord(str[i]);
    end;
  BKDRHash := res;
end;


function DEKHash(const str : string) : cardinal;
var
  i, res : cardinal;
begin
  res := length(str);
  for i := 1 to length(str) do
    begin
      res := ((res shr 5) xor (res shl 27)) xor ord(str[i]);
    end;
  DEKHash := res;
end;


function DJBHash(const str : string) : cardinal;
var
  i, res : cardinal;
begin
  res := 5381;
  for i := 1 to length(str) do
    begin
      res := ((res shl 5) + res) + ord(str[i]);
    end;
  DJBHash := res;
end;


procedure mk_hfp(fstr:fragstr);
var
  h1, h2, h3, h4 : cardinal;
begin
  h1 := (FNVHash(fstr) mod hfpsize) + 1;
  h2 := (DJBHash(fstr) mod hfpsize) + 1;
  if n_hfpbits > 2 then h3 := (BKDRHash(fstr) mod hfpsize) + 1;
  if n_hfpbits > 3 then h4 := (DEKHash(fstr) mod hfpsize) + 1;
  hfp[h1] := true;
  hfp[h2] := true;
  if n_hfpbits > 2 then hfp[h3] := true;
  if n_hfpbits > 3 then hfp[h4] := true;
end;


function assemble_frag(fp:ringpath_type;is_ring:boolean):fragstr;
var
  i, a1, a2, b : integer;
  pl : smallint;
  el1 : str2;
  fstr : string;
  bt, nbt : char;
  lfrag : frag_rec;
  newfragstr : fragstr;
  valid : boolean;
begin
  fstr := '';
  valid := true;
  fillchar(lfrag,sizeof(frag_rec),0);
  pl := path_length(fp);
  a1 := fp[1];
  el1 := atom^[a1].element;
  lfrag.size := 1;
  lfrag.a_atomicnum[1] := atomicnumber(el1);
  if (el1[2] = ' ') then delete(el1,2,1);
  fstr := fstr + el1;
  if pl > 1 then 
    begin
      lfrag.size := pl;
      lfrag.ring := is_ring;
      for i := 2 to pl do
        begin
          a2 := a1;
          a1 := fp[i];
          b := get_bond(a1,a2);
          el1 := atom^[a1].element;
          bt := bond^[b].btype;
          if bt = 'S' then nbt := '-';
          if bt = 'D' then nbt := '=';
          if (bt = 'A') or (bond^[b].arom = true) then nbt := '~';  // v0.4b  switched lines with the next one
          if bt = 'T' then nbt := '#';  // v0.4b  'T' remains 'T' even in aromatic rings
          if (pos(bt,'lsda')>0) then valid := false;  // reject fragments with wildcard bonds
          lfrag.a_atomicnum[i] := atomicnumber(el1);
          if lfrag.a_atomicnum[i] > 900 then valid := false;  // reject fragments with wildcard atoms
          lfrag.b_code[(i-1)] := nbt;
          if (el1[2] = ' ') then delete(el1,2,1);
          fstr := fstr + nbt + el1 ;
        end;
      if is_ring then   // add the last bond (between first and last atom)
        begin
          a1 := fp[pl];
          a2 := fp[1];
          b := get_bond(a1,a2);
          bt := bond^[b].btype;
          if bt = 'S' then nbt := '-';
          if bt = 'D' then nbt := '=';
          if bt = 'T' then nbt := '#';
          if (bt = 'A') or (bond^[b].arom = true) then nbt := '~';
          if (pos(bt,'lsda')>0) then valid := false;  // reject fragments with wildcard bonds
          lfrag.b_code[pl] := nbt;
          fstr := fstr + nbt;
        end;
    end;
  if valid then
    begin
      if (is_ring = false) then order_frag(lfrag);
      newfragstr := mk_fragstr(lfrag);
      assemble_frag := newfragstr;
    end else assemble_frag := '';
end;


procedure collect_frags(fp:ringpath_type);  // recursive procedure
var
  lfp : ringpath_type;
  nb : neighbor_rec;
  pl : smallint;
  n_nb : smallint;
  i, j, a_last : integer;
  thisfrag : fragstr;
  elem : str2;
begin
  pl := path_length(fp);
  fillchar(lfp,sizeof(ringpath_type),0);
  fillchar(nb,sizeof(neighbor_rec),0);
  lfp := fp;
  a_last := lfp[pl];
  if (pl > max_fragpath_length) then exit;
  if (atom^[a_last].heavy = false) then exit;
  // check if the last atom is a query atom and stop if yes
  elem := atom^[a_last].element;
  if (elem = 'Q ') or (elem = 'A ') or (elem = 'X ') then
    begin
      {$IFDEF debug}
      debugoutput('skipping query atoms in fragment search'); 
      {$ENDIF}
      lfp[pl] := 0;  // remove last atom (this could be the only one...)
      exit;
    end;
  // check if fragment forms a ring; if yes, discard it (rings are treated separately
  if (pl > 2) then
    begin
      for i := 1 to (pl-1) do
        begin
          if lfp[pl] = lfp[i] then exit;
        end;
    end;
  if (pl >= min_fragpath_length) and (pl <= max_fragpath_length) then
    begin
      //assemble the fragment and put it into the list
      thisfrag := assemble_frag(lfp,false);
      if (thisfrag <> '') then mk_hfp(thisfrag);
    end;
  if (pl = 1) then nb := get_neighbors(a_last) else nb := get_nextneighbors(a_last,lfp[(pl-1)]);
  // any other case: fragment path is not yet complete
  if (pl = 1) then n_nb := atom^[a_last].neighbor_count else n_nb := atom^[a_last].neighbor_count - 1;
  if (n_nb > 0) then
    begin
      for j := 1 to n_nb do
        begin
          lfp[(pl+1)] := nb[j];
          collect_frags(lfp);
        end;
    end;
end;  // collect_frags


procedure make_linear_frags;
var
  a : integer;
  lfp : ringpath_type;
begin
  if (n_atoms < 2) or (n_bonds < 1) then exit;
  fillchar(lfp,sizeof(ringpath_type),0);
  for a := 1 to n_atoms do
    begin
      lfp[1] := a;
      collect_frags(lfp);
    end;
end;


procedure rotate_ring(var rp:ringpath_type);
var
  i, a, rs : integer;
begin
  rs := path_length(rp);
  if (rs < 3) then exit;
  a := rp[rs];
  for i := rs downto 2 do rp[i] := rp[(i-1)];
  rp[1] := a;
end;


procedure reverse_ring(var rp:ringpath_type);
var
  i, a, rs : integer;
begin
  rs := path_length(rp);
  if rs < 3 then exit;
  for i := 1 to (rs div 2) do
    begin
      a := rp[i];
      rp[i] := rp[(rs-(i-1))];
      rp[(rs-(i-1))] := a;
    end;
end;


procedure make_ring_frags;
var
  i, j : integer;
  lfp : ringpath_type;
  h1, h2, h3, h4 : int64;
  fstr : fragstr;
  rs : integer;
begin
  if (n_rings < 1) then exit;
  for i := 1 to n_rings do
    begin
      rs := ringprop^[i].size;
      if (rs <= max_ringfragpath_length) and (ringprop^[i].envelope = false) then
        begin
          lfp := ring^[i];
          fstr := assemble_frag(lfp,true);
          if fstr = '' then exit;
          h1 := (FNVHash(fstr) mod hfpsize) + 1;
          h2 := (DJBHash(fstr) mod hfpsize) + 1;
          if n_hfpbits > 2 then h3 := (BKDRHash(fstr) mod hfpsize) + 1;
          if n_hfpbits > 3 then h4 := (DEKHash(fstr) mod hfpsize) + 1;
          //rotate the ring into every possible position
          for j := 1 to (rs-1) do 
            begin
              rotate_ring(lfp);
              fstr := assemble_frag(lfp,true);
              if fstr = '' then exit;
              h1 := h1 + ((FNVHash(fstr) mod hfpsize) + 1);
              h2 := h2 + ((DJBHash(fstr) mod hfpsize) + 1);
              if n_hfpbits > 2 then h3 := h3 + ((BKDRHash(fstr) mod hfpsize) + 1);
              if n_hfpbits > 3 then h4 := h4 + ((DEKHash(fstr) mod hfpsize) + 1);
            end;

          //now reverse the ring and repeat the whole procedure
          lfp := ring^[i];
          reverse_ring(lfp);
          fstr := assemble_frag(lfp,true);
          if fstr = '' then exit;
          h1 := h1 + ((FNVHash(fstr) mod hfpsize) + 1);
          h2 := h2 + ((DJBHash(fstr) mod hfpsize) + 1);
          if n_hfpbits > 2 then h3 := h3 + ((BKDRHash(fstr) mod hfpsize) + 1);
          if n_hfpbits > 3 then h4 := h4 + ((DEKHash(fstr) mod hfpsize) + 1);
          for j := 1 to (rs-1) do 
            begin
              rotate_ring(lfp);
              fstr := assemble_frag(lfp,true);
              if fstr = '' then exit;
              h1 := h1 + ((FNVHash(fstr) mod hfpsize) + 1);
              h2 := h2 + ((DJBHash(fstr) mod hfpsize) + 1);
              if n_hfpbits > 2 then h3 := h3 + ((BKDRHash(fstr) mod hfpsize) + 1);
              if n_hfpbits > 3 then h4 := h4 + ((DEKHash(fstr) mod hfpsize) + 1);
            end;

          //accumulated hash values have to be broken down again
          //writeln(' sum: ',h1,'  ',h2);
          h1 := (h1 mod hfpsize) + 1;
          h2 := (h2 mod hfpsize) + 1;
          if n_hfpbits > 2 then h3 := (h3 mod hfpsize) + 1;
          if n_hfpbits > 3 then h4 := (h4 mod hfpsize) + 1;
          
          //writeln('==cumulative ring hash for ring ',i,' ========> ',h1,', ',h2);

          hfp[h1] := true;
          hfp[h2] := true;
          if n_hfpbits > 2 then hfp[h3] := true;
          if n_hfpbits > 3 then hfp[h4] := true;
        end;
    end;
end;

{$IFDEF extended_hfp}
// v0.4e; added extra hfp bits for branched carbons and nitrogens

procedure reorder_frag(var frag:frag_rec);
// order branched fragments by atomic number + bond type
var
  i, j, k : integer;
  score1, score2 : longint;
  lfrag : frag_rec;
  newpos : integer;
begin
  fillchar(lfrag,sizeof(frag_rec),0);
  lfrag.size := 0;
  lfrag.ring := frag.ring;
  for i := 1 to frag.size do
    begin
      score1 := frag.a_atomicnum[i] * 10;
      case frag.b_code[i] of
        '-' : inc(score1,1);
        '=' : inc(score1,2);
        '#' : inc(score1,3);
        '~' : inc(score1,4);
      end;
      inc(lfrag.size);
      if (lfrag.size = 1) then 
        begin
          lfrag.a_atomicnum[1] := frag.a_atomicnum[1];
          lfrag.b_code[1] := frag.b_code[1];          
        end else
        begin
          // sort any additional entry
          newpos := 1;
          for j := 1 to (lfrag.size -1) do
            begin
              score2 := lfrag.a_atomicnum[j] * 10;
              case lfrag.b_code[j] of
                '-' : inc(score2,1);
                '=' : inc(score2,2);
                '#' : inc(score2,3);
                '~' : inc(score2,4);
              end;
              if (score2 > score1) then inc(newpos);
            end;
          if (newpos < lfrag.size) then
            begin
              for k := lfrag.size downto (newpos+1) do
                begin
                  lfrag.a_atomicnum[k] := lfrag.a_atomicnum[k-1];
                  lfrag.b_code[k] := lfrag.b_code[k-1];
                end;
            end;
          lfrag.a_atomicnum[newpos] := frag.a_atomicnum[i];
          lfrag.b_code[newpos] := frag.b_code[i];
        end;
    end;  // for i ...
  // now copy the ordered elements back to frag
  for i := 1 to frag.size do
    begin
      frag.a_atomicnum[i] := lfrag.a_atomicnum[i];
      frag.b_code[i] := lfrag.b_code[i];
    end;
end;

function assemble_brfragstr(center_atomicnum:integer;frag:frag_rec):fragstr;
var
  i, code : integer;
  fstr, tmpstr : string;
begin
  fstr := '';
  fstr := inttostr(center_atomicnum) + '(';
  for i := 1 to frag.size do
    begin
      if (i > 1) then fstr := fstr + ',';
      fstr := fstr + frag.b_code[i];
      fstr := fstr + inttostr(frag.a_atomicnum[i]);
    end;
  fstr := fstr + ')';
  assemble_brfragstr := fstr;
end;

procedure make_branched_frags;
var
  i, j, k, a1, a2, b, can : integer;
  nb : neighbor_rec;
  n_nb : smallint;
  thisfrag : fragstr;
  elem, el2 : str2;
  valid : boolean;
  lfrag, lfrag3 : frag_rec;
  bt, nbt : char;
  brfragstr : fragstr;
begin
  if (n_atoms > 3) then
    begin
      for i := 1 to n_atoms do
        begin
          a1 := i;
          elem := atom^[i].element;        
          n_nb := atom^[i].neighbor_count; 
          if (((elem = 'C ') or (elem = 'N ')) and ((n_nb = 4) or (n_nb = 3))) then
            begin
              valid := true;
              fillchar(nb,sizeof(neighbor_rec),0);
              lfrag.size := n_nb;
              fillchar(lfrag.a_atomicnum,max_ringfragpath_length,0);
              fillchar(lfrag.b_code,max_ringfragpath_length,' ');
              lfrag.ring := false;
              nb := get_neighbors(i);
              can := atomicnumber(elem);
              for j := 1 to n_nb do
                begin
                  a2 := nb[j];
                  b := get_bond(a1,a2);
                  bt := ' ';
                  nbt := ' ';
                  if (b > 0) then bt := bond^[b].btype else valid := false;
                  el2 := atom^[a2].element;
                  if bt = 'S' then nbt := '-';
                  if bt = 'D' then nbt := '=';
                  if (bt = 'A') or (bond^[b].arom = true) then nbt := '~';
                  if bt = 'T' then nbt := '#';
                  if (pos(bt,'lsda')>0) then valid := false;  // reject fragments with wildcard bonds
                  lfrag.a_atomicnum[j] := atomicnumber(el2);
                  if lfrag.a_atomicnum[j] > 900 then valid := false;  // reject fragments with wildcard atoms
                  if ((el2 = 'A ') or (el2 = 'Q ') or (el2 = 'X ')) then valid := false;  // reject fragments with wildcard atoms
                  lfrag.b_code[j] := nbt;
                end;  // for j...
              if valid then
                begin
                  //brfragstr := assemble_brfragstr(can,lfrag);
                  //writeln('  branched atom:',i,'  ',brfragstr);
                  reorder_frag(lfrag);
                  brfragstr := assemble_brfragstr(can,lfrag);
                  //writeln('        primary fragment after reordering: ',brfragstr);
                  mk_hfp(brfragstr);
                  // if this is a quaternary C or N, make all possible tertiary subgraphs
                  if (n_nb = 4) then
                    begin
                      for j := 1 to 4 do  // outer loop determines which substituent NOT to include
                        begin
                          fillchar(lfrag3,sizeof(frag_rec),0);
                          lfrag3.size := 0;
                          lfrag3.ring := false;
                          for k := 1 to 4 do  // inner loop copies the remaining substituents from lfrag to lfrag3
                            begin
                              if (j <> k) then
                                begin
                                  inc(lfrag3.size);
                                  lfrag3.a_atomicnum[lfrag3.size] := lfrag.a_atomicnum[k];
                                  lfrag3.b_code[lfrag3.size] := lfrag.b_code[k];
                                end;
                            end;
                          // and now the hash

                          reorder_frag(lfrag3);
                          brfragstr := assemble_brfragstr(can,lfrag3);
                          //writeln('        sub-branch after reordering: ',brfragstr);
                          mk_hfp(brfragstr);
                        end;
                    end;
                end;  // if valid
            end;  // elem = C ...
        end;  // for i...
    end;
end;

{$ENDIF}


procedure write_hfp;
var
  i, n1 : integer;
begin
  n1 := 0;
  for i := 1 to hfpsize do
    begin
      if hfp[i] then 
        begin
          write('1');
          inc(n1);
        end else write('0');
    end;
  write(';',n1);
  if opt_verbose then write('  darkness (%): ',(100*n1/hfpsize):1:1);
  writeln;
end;


procedure write_hfp_dec;
const
  bsize = 32;
var
  i, j, n1 : integer;
  hfpincrement : int64;
  hfpdecimal : int64;
begin
  n1 := 0;
  for i := 0 to ((hfpsize div bsize) - 1) do
    begin
      hfpdecimal := 0;
      for j := 1 to bsize do
        begin
          hfpincrement := 1;
          if hfp[((bsize*i)+j)] then 
            begin
              inc(n1);
              hfpincrement := hfpincrement shl (j-1);
              hfpdecimal := hfpdecimal + hfpincrement;
            end;
        end;
      if (i > 0) then write(',');
      write(hfpdecimal);
    end;
  write(';',n1);
  writeln;
end;


procedure make_hashed_fp;
begin
  init_pt;
  fillchar(hfp,sizeof(hfp),false);
  // get the fragments
  make_linear_frags;
  make_ring_frags;
  {$IFDEF extended_hfp}
  make_branched_frags;
  {$ENDIF}
  // now write the output
  if hfpformat = fpf_boolean then write_hfp;
  if hfpformat = fpf_decimal then write_hfp_dec; 
end;


procedure roll_fwd(var xnb:neighbor_rec;steps:integer);  // v0.4a
var
  i,t,n : integer;
  xnb_size : integer;
begin
  if (steps < 1) then exit;
  xnb_size := 0;
  for i := 1 to max_neighbors do if (xnb[i] > 0) then xnb_size := i;
  for n := 1 to steps do
    begin
      t := xnb[1];
      for i := 2 to xnb_size do xnb[(i-1)] := xnb[i];
      xnb[xnb_size] := t;
    end;
end;

function count_gmm_matches:integer;  // v0.4a
var
  i, j, n : integer;
begin
  n := 0;
  for i := 1 to ndl_n_atoms do
    begin
      for j := 1 to n_atoms do if gmm^[i,j] then inc(n);
    end;
  count_gmm_matches := n;
end;

function gmm_collision:boolean;   // v0.4b
// scans global match matrix for collisions: i.e., matches of a single haystack
//atom with (at least) two different needle atoms which do not have an alternative
// match (e.g., 4-methylheptane would match 1,4-dimethylcyclohexane
var
  i, j, k, l : integer;
  n_hits, hit_pos : integer;
  n_alt, n_coll, dup_pos : integer;
  res : boolean;
begin
  if not (use_gmm and valid_gmm) then
    begin
      gmm_collision := false;
      exit;
    end;
  n_coll := 0;
  res := false;
  for i := 1 to ndl_n_atoms do
    begin
      n_hits  := 0;
      hit_pos := 0;
      for j := 1 to n_atoms do
        begin
          if gmm^[i,j] then
            begin
              inc(n_hits);
              hit_pos := j;
            end;
        end;
      if (n_hits = 1) then
        begin
          dup_pos := 0;
          for k := 1 to ndl_n_atoms do
            begin
              if (k <> i) and gmm^[k,hit_pos] then 
                begin
                  dup_pos := k;;
                  n_alt := 0;  // number of alternative matches
                  for l := 1 to n_atoms do
                    begin
                      if (l <> hit_pos) and gmm^[dup_pos,l] then inc(n_alt);
                    end;
                  if n_alt < 1 then inc(n_coll);
                end;
            end;  
        end;
        
    end;
  if (n_coll > 0) then res := true;
  gmm_collision := res;
end;


procedure cleanup_gmm;  // v0.4a
var
  i, j, k, l : integer;
  old_n_matches, new_n_matches : integer;
  n_hits, hit_pos : integer;
  ndl_nb : neighbor_rec;
  hst_nb : neighbor_rec;
  ndl_n_nb, hst_n_nb : integer;
  na, ha : integer;
  found_neighbor : boolean;
  found_all_neighbors : boolean;
begin
  new_n_matches := count_gmm_matches;
  repeat
    old_n_matches := new_n_matches;
    // first step: check for unique matches and remove other substituents from here
    for i := 1 to ndl_n_atoms do
      begin
        n_hits := 0;
        for j := 1 to n_atoms do
          begin
            if gmm^[i,j] then
              begin
                inc(n_hits);
                hit_pos := j;
              end;
          end;
        if (n_hits = 1) then
          begin
            for k := 1 to ndl_n_atoms do
              begin
                if (k <> i) then gmm^[k,hit_pos] := false;
              end;  
          end;
      end;
    // second step: check for dangling neighbor atoms
    for i := 1 to ndl_n_atoms do
      begin
        fillchar(ndl_nb,sizeof(neighbor_rec),0);  // v0.4e: moved this here from inner looüp
        ndl_nb := get_ndl_neighbors(i);
        ndl_n_nb := ndl_atom^[i].neighbor_count;
        for j := 1 to n_atoms do
          begin
            if gmm^[i,j] then
              begin
                fillchar(hst_nb,sizeof(neighbor_rec),0);
                hst_nb := get_neighbors(j);
                hst_n_nb := atom^[j].neighbor_count;
                if (ndl_n_nb > 0) then
                  begin
                    found_all_neighbors := true;
                    for k := 1 to ndl_n_nb do
                      begin
                        found_neighbor := false;
                        for l := 1 to hst_n_nb do
                          begin
                            na := ndl_nb[k];
                            ha := hst_nb[l];
                            if gmm^[na,ha] then 
                              begin 
                                found_neighbor := true;
                              end;
                          end;
                        if (found_neighbor = false) then found_all_neighbors := false;
                      end;
                    if (found_all_neighbors = false) then gmm^[i,j] := false;
                  end;
              end;          
          end;
      end;
    new_n_matches := count_gmm_matches;
  until new_n_matches = old_n_matches;
end;

procedure copy_gmm_to_total;  // v0.4a
var
  i, j : integer;
begin
  for i := 1 to ndl_n_atoms do
    for j := 1 to n_atoms do
      if gmm^[i,j] then gmm_total^[i,j] := true;
end;

(*
procedure write_matches;  // v0.4a
var
  i,j,k,nx,npos : integer;
  mnb : neighbor_rec;
begin
  write(' ');
  for i := 1 to ndl_n_atoms do
    begin
      nx := 0;
      fillchar(mnb,sizeof(neighbor_rec),0);
      for j := 1 to n_atoms do
        begin
          if gmm_total^[i,j] then
            begin
              inc(nx);
              mnb[nx] := j;
            end;
        end;
      if (nx > 0) then
        begin
          if (i > 1) then write(';');
          write(i,'=');
          if (nx > 1) then
            begin
              npos := 0;
              for k := 1 to i do
                begin
                  if gmm_total^[k,(mnb[1])] then inc(npos);
                  //if gmm_total^[k,(mnb[1])] then write('-X-') else write('-o-');
                end;
              if (npos > 1) then roll_fwd(mnb,(npos-1));
              if opt_matchnum then
                begin
                  for k := 1 to nx do
                    begin
                      if (k > 1) then write(',');
                      write(mnb[k]);
                    end;
                  end;
              if opt_matchnum1 then
                begin
                  write(mnb[1]);
                end;
            end else write(mnb[1]);
        end;
    end;
end;
*)

procedure write_matches;  // v0.4a, v0.4b
var
  i,j,k,nx,npos : integer;
  mnb : neighbor_rec;
begin
  if (length(match_string) > 0) then match_string := match_string + '.';
  for i := 1 to ndl_n_atoms do
    begin
      nx := 0;
      fillchar(mnb,sizeof(neighbor_rec),0);
      for j := 1 to n_atoms do
        begin
          if gmm_total^[i,j] then
            begin
              inc(nx);
              mnb[nx] := j;
            end;
        end;
      if (nx > 0) then
        begin
          if (i > 1) then match_string := match_string +';';
          match_string := match_string + inttostr(i) + '=';
          if (nx > 1) then
            begin
              npos := 0;
              for k := 1 to i do
                begin
                  if gmm_total^[k,(mnb[1])] then inc(npos);
                  //if gmm_total^[k,(mnb[1])] then write('-X-') else write('-o-');
                end;
              if (npos > 1) then roll_fwd(mnb,(npos-1));
              if opt_matchnum then
                begin
                  for k := 1 to nx do
                    begin
                      if (k > 1) then match_string := match_string +',';
                      match_string := match_string + inttostr(mnb[k]);
                    end;
                  end;
              if opt_matchnum1 then
                begin
                  match_string := match_string + inttostr(mnb[1]);
                end;
            end else match_string := match_string + inttostr(mnb[1]);
        end;
    end;
end;


begin  // main routine
  progname := extractfilename(paramstr(0));
  if pos('MATCHMOL',upcase(progname))>0 then progmode := pmMatchMol else
    begin
      if pos('CHECKMOL',upcase(progname))>0 then progmode := pmCheckMol else
        begin
          writeln('THOU SHALLST NOT RENAME ME!');
          halt(9);
        end;
    end;
  if (paramcount = 0) then
    begin
      show_usage;
      halt(1);
    end;
  init_globals;
  init_molstat(molstat);
  parse_args;
  if ringsearch_mode = rs_sar then max_vringsize := max_ringsize else
                                   max_vringsize := ssr_vringsize;  // v0.3n (was: 10)
  //if opt_verbose then writeln(progname+' v',version,'  N. Haider 2003-2007');
  if progmode = pmMatchMol then
    begin
      left_trim(ndl_molfilename);
      left_trim(molfilename);
      if (((molfilename = '') or (ndl_molfilename = '') or (paramcount < 2))
        and (not opt_stdin)) then
        begin
          show_usage;
          halt(2);      // new in v0.2k
        end;
      if ((not fileexists(ndl_molfilename))  and (not opt_stdin)) then
        begin
          if ((length(ndl_molfilename) > 1) and (ndl_molfilename[1] = '-')) then
            show_usage else writeln('file ',ndl_molfilename,' not found!');  // new in v0.2k
          halt(2);
        end;
    end;
  if ((not fileexists(molfilename))  and (not opt_stdin)) then
    begin
      if ((length(molfilename) > 1) and (molfilename[1] = '-')) then
        begin
          //writeln('  you entered this molfilename: ',molfilename);
          show_usage;
        end else writeln('file ',molfilename,' not found!');  // new in v0.2k
      halt(2);
    end;
  // read the first molecule and process it; if we are in "matchmol" mode,
  // this is the "needle"
  if progmode = pmMatchMol then readinputfile(ndl_molfilename) else
                                readinputfile(molfilename);
  li := 1;   // initialize line pointer for input buffer
  filetype := get_filetype(ndl_molfilename);
  if filetype = 'unknown' then
    begin
      writeln('unknown query file format!');
      if opt_verbose then
        begin
          writeln('===========================================');
          for i := 1 to molbufindex do writeln(molbuf^[i]);
        end;
      halt(3);
    end;
  mol_OK := true;           // added in v0.2i
  if filetype = 'alchemy' then read_molfile(ndl_molfilename);
  if filetype = 'sybyl'   then read_mol2file(ndl_molfilename);
  if filetype = 'mdl'     then read_MDLmolfile(ndl_molfilename);
  count_neighbors;
  if (not mol_OK) or (n_atoms < 1) then  // v0.3g; check if this is a valid query structure
    begin
      writeln('invalid molecule');
      halt(3);
    end;
  if (not found_arominfo) or (progmode = pmCheckMol) then  // added in v0.2b/0.2c
    begin
      {$IFDEF debug}
      if (not found_arominfo) then debugoutput('no aromaticity information found - checking myself...')
        else debugoutput('performing full aromaticity check');   // new in v0.3d
      {$ENDIF}
      chk_ringbonds;
      if ringsearch_mode = rs_ssr then remove_redundant_rings;
      if n_rings = max_rings then
        begin
          if opt_verbose then writeln('warning: max. number of rings exceeded, reverting to SSR search');
          ringsearch_mode := rs_ssr;
          auto_ssr := true;  // v0.3n
          clear_rings;
          max_vringsize := ssr_vringsize;   // v0.3n (was: 10)
          chk_ringbonds;
          remove_redundant_rings;
        end;
      update_ringcount;
      // new in v0.3k: if output is a molfile, leave the original
      // representation of N-oxides, S-oxides, nitro groups, etc.
      // unchanged (ionic or non-ionic), in any other case make covalent bonds
      if (not opt_xmdlout) then normalize_ionic_bonds;  // v0.3k
      update_atypes;
      update_Htotal;    // added in v0.3
      chk_arom;
      if (ringsearch_mode = rs_ssr) then   // new in v0.3
        begin
          repeat
            prev_n_ar := count_aromatic_rings;
            chk_arom;
            n_ar := count_aromatic_rings;
          until ((prev_n_ar - n_ar) = 0);
        end;
    end else
    begin
      {$IFDEF debug}
      debugoutput('found aromaticity information in input file');
      {$ENDIF}
      if (not opt_xmdlout) then normalize_ionic_bonds;  // v0.3k
      update_atypes;  // added in v0.2f
      update_Htotal;  // end v0.2b snippet
    end;
  if progmode = pmCheckMol then
    begin
      if opt_verbose then write_mol;
      get_molstat;
      if (opt_molstat or opt_hfp) then
        begin
          if opt_molstat then
            begin
              if opt_molstat_X then write_molstat_X else write_molstat;
            end;
        end else
        begin
          if found_querymol then
            begin
              writeln('input structure contains query atom or query bond!');
              halt(1);
            end;
        end;
      if opt_none then opt_text := true;
      if (opt_text or opt_text_de or opt_code or opt_bin or opt_bitstring or opt_pos) then  // v0.5
        begin  
          chk_functionalgroups;
          if opt_text      then write_fg_text(lang_en);
          if opt_text_de   then write_fg_text(lang_de);
          if opt_code      then write_fg_code;
          if opt_bin       then write_fg_binary;
          if opt_bitstring then write_fg_bitstring;
          if opt_pos       then write_fg_pos;  // v0.5
        end;
      if opt_xmdlout   then 
        begin
          init_pt;
          calc_mf_mw;  // v0.4d
          write_MDLmolfile;
        end;
      if opt_hfp then make_hashed_fp;
      //if opt_verbose   then write_mol;
      zap_molecule;
    end else
    begin  // progmode = pmMatchMol
      // v0.4a
      if use_gmm then
        begin
          try
            getmem(gmm,sizeof(global_matchmatrix));
            getmem(gmm_total,sizeof(global_matchmatrix));
          except
            on e:Eoutofmemory do
              begin
                writeln('Not enough memory for global match matrix');
                halt(4);
              end;
          end;
        end;
      // now transfer all data to the "needle" set of variables, except for "fingerprint" mode
      if (not opt_fp) then   // v0.3m
        begin
          copy_mol_to_needle;
          chk_wildcard_rings;  // v0.3p
          set_ndl_atom_tags;   // v0.3o
          if opt_verbose then write_needle_mol;
          if (ndl_n_atoms > max_ndl_gmmsize) then valid_gmm := false else valid_gmm := true;  // v0.4a
          if rs_strict then      // v0.3j
            begin
              ndl_ref_atom := find_ndl_ref_atom_cv;
            end
          else
            ndl_ref_atom := find_ndl_ref_atom;
        end else 
        begin
          copy_mol_to_tmp;    // v0.3m
          if opt_verbose then writeln('1st molecule stored in buffer: ',tmp_molname);
          valid_gmm := false;  // v0.4a global match matrix is not needed in fingerprint mode
        end;
      // next, read the "haystack" file and process it
      li := 1;
      mol_count := 0;
      fpdecimal := 0;  // v0.3m
      fpindex   := 0;  // v0.3m
      repeat
        begin
          // new in v0.3i: reset ringsearch_mode to its initial value
          // for each new molecule
          ringsearch_mode := opt_rs;
          if ringsearch_mode = rs_sar then max_vringsize := max_ringsize else
                                           max_vringsize := ssr_vringsize;  // v0.3n (was: 10)
          readinputfile(molfilename);
          li := 1;
          filetype := get_filetype(molfilename);
          if filetype <> 'unknown' then
            begin
              found_arominfo := false;  // added in v0.2b
              mol_OK := true;           // added in v0.2i
              if filetype = 'alchemy' then read_molfile(molfilename);
              if filetype = 'sybyl'   then read_mol2file(molfilename);
              if filetype = 'mdl'     then read_MDLmolfile(molfilename);
              inc(mol_count);
              inc(fpindex);
              count_neighbors;
              //if (not mol_OK) or (n_atoms < 1) then writeln(mol_count,':no valid structure found') else
              if (not mol_OK) or (n_atoms < 1) and (not (opt_fp and (fpformat = fpf_decimal))) then 
                writeln(mol_count,':F') else    // v0.3l
                begin
                  if opt_exact and (not ndl_querymol) and   // v0.3p
                     ((n_Ctot <> ndl_n_Ctot) or (n_Otot <> ndl_n_Otot) or
                     (n_Ntot <> ndl_n_Ntot)) then       // new in v0.3g
                    begin
                      if (not opt_molout) and (not (opt_fp and (fpformat = fpf_decimal))) then writeln(mol_count,':F');
                    end else
                    begin 
                      if ((not found_arominfo) or (opt_strict and tmfmismatch)) then  // added in v0.3m
                        begin
                          {$IFDEF debug}
                          debugoutput('no aromaticity information found (or tweak mismatch) - checking myself...');
                          {$ENDIF}
                          chk_ringbonds;
                          if ringsearch_mode = rs_ssr then remove_redundant_rings;
                          if n_rings = max_rings then
                            begin
                              if opt_verbose then writeln('Warning: max. number of rings reached, reverting to SSR search');
                              ringsearch_mode := rs_ssr;
                              clear_rings;
                              max_vringsize := ssr_vringsize;  // v0.3n (was: 10)
                              chk_ringbonds;
                              remove_redundant_rings;
                            end;
                          update_ringcount;
                          update_atypes;
                          update_Htotal;   // added in v0.3
                          chk_arom;
                          if (ringsearch_mode = rs_ssr) then   // new in v0.3
                            begin
                              repeat
                                prev_n_ar := count_aromatic_rings;
                                chk_arom;
                                n_ar := count_aromatic_rings;
                              until ((prev_n_ar - n_ar) = 0);
                            end;
                        end else 
                        begin
                          {$IFDEF debug}
                          debugoutput('found aromaticity information in input file');
                          {$ENDIF}
                          if opt_strict then update_atypes;  // added in v0.2f
                          update_Htotal;
                        end;
                      init_molstat(ndl_molstat);
                      if normalize_ionic_bonds then update_atypes;   // new in v0.3k, modified in v0.3m
                      // if in "fingerprint mode", exchange needle and haystack
                      if (opt_fp) then   // v0.3m
                        begin
                          zap_needle; 
                          copy_mol_to_needle;
                          chk_wildcard_rings;  // v0.3p
                          zap_molecule;
                          copy_tmp_to_mol;
                          if opt_verbose then write_needle_mol;
                          if rs_strict then      // v0.3j
                            ndl_ref_atom := find_ndl_ref_atom_cv
                          else
                            ndl_ref_atom := find_ndl_ref_atom;
                          if opt_verbose then write_mol;
                        end;    // v0.3m
                      // now that we have both molecules, perform the comparison
                      // v0.3o: takes care of disconnected fragment...
                      update_atypes_quick;  // v0.4b
                      if opt_verbose and (not opt_fp) then write_mol;  // v0.4b
                      clear_atom_tags;
                      overall_match := false;  // v0.4a
                      fillchar(tmp_tag,sizeof(tmp_tag),false);   // v0.4a
                      found_untagged := false;
                      if (use_gmm and valid_gmm) then fillchar(gmm_total^,sizeof(global_matchmatrix),false);  // v0.4a
                      n_matches := 0;  // v0.4b
                      repeat          // this repeat...until is new in v0.4a
                        begin
                          set_ndl_atom_tags;
                          matchsummary := true;
                          if (use_gmm and valid_gmm) then clear_gmm;  // v0.4a
                          perform_match;
                          matchsummary := matchresult;
                          if (count_tagged_ndl_heavyatoms > 0) and (matchsummary = true) then
                            begin
                              repeat
                                begin
                                  if rs_strict then
                                    ndl_ref_atom := find_ndl_ref_atom_cv
                                  else
                                    ndl_ref_atom := find_ndl_ref_atom;
                                  perform_match;
                                  if (matchresult = false) then matchsummary := false;
                                end;
                              until (count_tagged_ndl_heavyatoms = 0) or (matchsummary = false);  
                            end;
                          // end of disconnected-fragment matching (v0.3o)
                          if matchsummary then    // v0.4b
                            begin
                              overall_match := true;
                              inc(n_matches);
                              if (use_gmm and valid_gmm) then
                                begin
                                  cleanup_gmm;
                                  clear_gmm_total;  // v0.4b
                                  copy_gmm_to_total;
                                  write_matches;   // v0.4b
                                end;
                            end;
                        end;
                      until ((use_gmm = false) or (valid_gmm = false) or (n_matches >= max_n_matches) or  // v0.4b
                            (opt_matchnum1 = true) or (found_untagged = false));  // end v0.4a
                      if overall_match = true then  // v0.3o, v0.4a
                        begin
                          if opt_molout then
                            begin
                              for i := 1 to molbufindex do writeln(molbuf^[i]);
                              if opt_verbose and use_gmm then write_gmm;   // v0.4a
                            end else 
                            begin
                              if (not opt_fp) then
                                begin
                                  write(inttostr(mol_count),':T');   // v0.4a
                                  if (use_gmm and valid_gmm and (opt_matchnum or opt_matchnum1)) then 
                                    begin
                                      //write_matches;  // v0.4a
                                      write(' ');           // v0.4b
                                      write(match_string);  // v0.4b
                                      match_string := '';   // v0.4b
                                    end;
                                  writeln(' '+sep_label);       // v0.4a, v0.4b
                                  sep_label := '';          // v0.4b
                                end else
                                begin
                                  //if (ndl_n_heavyatoms = n_heavyatoms) and
                                  //   (ndl_n_heavybonds = n_heavybonds) then 
                                  //   fp_exacthit := true else fp_exacthit := false;
                                  if (ndl_n_atoms = n_atoms) and   // v0.4b; also explicit H must match!
                                     (ndl_n_bonds = n_bonds) then 
                                     fp_exacthit := true else fp_exacthit := false;
                                  if fp_exacthit then fp_exactblock := true;
                                  if (fpformat = fpf_boolean) then
                                    begin
                                      if fp_exacthit then
                                        writeln(inttostr(mol_count),':TX') 
                                      else
                                        writeln(inttostr(mol_count),':T');
                                    end;
                                  if (fpformat = fpf_decimal) then
                                    begin
                                      fpincrement := 1;
                                      for i := 1 to fpindex do fpincrement := fpincrement shl 1;
                                      fpdecimal := fpdecimal + fpincrement;                        
                                    end;                                     
                                end;
                            end;
                        end else
                        begin
                          if not (opt_molout or (opt_fp and (fpformat = fpf_decimal))) then writeln(inttostr(mol_count),':F');
                        end;
                      if (opt_fp and (fpformat = fpf_decimal)) and (fpindex = fp_blocksize) then
                        begin
                          if fp_exactblock then fpdecimal := fpdecimal + 1;
                          writeln(fpdecimal);
                          fpindex   := 0;
                          fpdecimal := 0;
                          fp_exactblock := false;
                        end;
                      zap_molecule;
                      molbufindex := 0;
                    end;
                end;  // mol_OK
            end else writeln(mol_count,':unknown file format'); // if filetype <> 'unknown'
        end;
      until (mol_in_queue = false);
      if (opt_fp and (fpformat = fpf_decimal)) and (fpindex > 0) then
        begin
          if fp_exactblock then fpdecimal := fpdecimal + 1;
          writeln(fpdecimal);
        end;
      zap_needle;
      if rfile_is_open then close(rfile);  // new in v0.2g
      if use_gmm then   // v0.4a
        begin
          try
            freemem(gmm,sizeof(global_matchmatrix));
            freemem(gmm_total,sizeof(global_matchmatrix));
          except
          end;   
        end;   
    end;
end.

File organization and naming rules
==================================

The files are organized into three main folders: Gr_Kerr, MGr_dω_k, and HDG_QNMs.


1. Gr_Kerr
-------------

Folder structure:
  Gr_Kerr/
    l2/
    l3/
    l4/

Here, l2, l3, and l4 correspond to l = 2, 3, and 4.

File naming rule:
  l{l}m{m}n{n}or{WKB order}_{re/im}.txt

Examples:
  l2m0n0or1_im.txt
    means l = 2, m = 0, n = 0, WKB order = 1 and imaginary part
  l2mm2n2or2_re.txt
    means l = 2, m = -2, n = 2, WKB order = 2 and real part
    Here, "mm2" means m = -2.

File contents:
  Each row contains values of:
  a     Re/Im[ω] 
  where a is the spin parameter, Re/Im[ω]  is real/imaginary part of QNM frequency


2. MGr_dω_k
------------------

Folder structure:
  MGr_dω_k/
    l2n0/

Here, l2n0 corresponds to l = 2 and n = 0

File naming rule:
  l{l}m{m}n{n}or{WKB order}_{re/im}.txt

Example:
  l2m2n0or3_im.txt
    means l = 2, m = 2, n = 0, WKB order = 3  and imaginary part

File contents:
  Each row contains values of:

    a      Re/Im[ω]      Re/Im[ω_k=4]      Re/Im[ω_k=3]      Re/Im[ω_k=2]      Re/Im[ω_k=1]      Re/Im[ω_k=0]      Re/Im[ω_k=-1]


3. HDG_QNMs
-----------

Folder structure:
  HDG_QNMs/
    l2m1n0/
    l2m2n0/

Here, for example, l2m1n0 corresponds to l = 2, m = 1, n = 0.

File naming rule:
  l{l}m{m}n{n}or{WKB order}Lam{lambda}_{even/odd}_{minus/plus}__{re/im}.txt

Example:
  l2m1n0or1Lam0.1_ev_minus_re.txt
    means l = 2, m = 1, n = 0, WKB order = 1, λ = 0.1,even, minus, real part

File contents:
  Each row contains values of:

   a      Re/Im[ω_HDG]     Re/Im[δω]

  where a is the spin parameter, Re/Im[ω_HDG] is the real/imaginary part of QNM frequency in Higher-Derivative Gravity, and Re/Im[δω] is the real/imaginary part of deviation from the corresponding Kerr value.

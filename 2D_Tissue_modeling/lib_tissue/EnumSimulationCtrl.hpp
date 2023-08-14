/* ENUMSIMULATIONCTRL_HPP
*	enum type incorporating Atrial regional cell types,
*	IKur mutation types
*	and AF types.

*   Haibo NI qiangzi.ni@gmail.com
*   Date: 21 January 2015 (Wednesday)
    15:04
*/


#ifndef ENUMSIMULATIONCTRL_HPP
#define ENUMSIMULATIONCTRL_HPP

#include <iostream>
#include <cstdlib>

enum TypeCell {RA = 0, PM, CT, RAA, AS, AVR, BB, LA, LAA, PV, SAN_C, SAN_P, AVN, LVEPI, LVMCELL, LVENDO, RVEPI, RVMCELL, RVENDO, PFI, PFMB, PF};
enum IKurMuationType {WT = 0, D322H, E48G, A305T, Y155C, D469E, P488S};
enum AFType {NONE = 0, AF1, AF2, AF3, AF4};




typedef TypeCell AtriaCellType;
TypeCell GetCellTypeFromLabel(int tissue);

#endif
// print R messages or warnings from C++
//
// Copyright (C) 2020 Karl Browman 
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// ###################################################################
// DISCLAIMER ON COPYRIGHT
// ###################################################################
//
// The above copyright notice was added by Robert Vogel to clearly
// distinguish the author's work, and license for which it is published, 
// from my own.  The source was copied from GitHub repository 
//
// https://github.com/rqtl/qtl2
//
// The DESCRIPTION file states that the author and current (April 2026)
// maintainer is Karl Broman, as such I have attributed the copyright
// to him despite not definitively seeing a statement of copyright.
// Moreover, I specified the year 2020 as, according to GitHub, the file
// was last modified 6 years ago.
//
// ###################################################################
// END DISCLAIMER
// ###################################################################
//
#ifndef R_MESSAGE_H
#define R_MESSAGE_H

#include <string>
#include <Rcpp.h>
#include <R_ext/Error.h>

#define RQTL2_NODEBUG 1 // ignore debugging code

void r_message(std::string text);
void r_warning(std::string text);

// Following based on code from Luke Miratrix, http://bit.ly/rcpp_assert
#ifdef RQTL2_NODEBUG
#define r_assert(EX)
#else
#define r_assert(EX) (void)((EX) || (__r_assert (#EX, __FILE__, __LINE__),0))
#endif

// r_enforce is like r_assert but always works
#define r_enforce(EX) (void)((EX) || (__r_assert (#EX, __FILE__, __LINE__),0))

#endif // R_MESSAGE_H

/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * ool_version.h
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * Sergio Drumond Ventura
 * Luis Alberto D'Afonseca
 * since: Feb, 28, 2004
 *
 * $Id: ool_version.h.in,v 1.1.1.1 2004/05/18 17:25:42 akiles Exp $
 *----------------------------------------------------------------------------*/

#ifndef __OOL_VERSION_H__
#define __OOL_VERSION_H__

#include <gsl/gsl_types.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS
  
#define OOL_VERSION "0.2.0"
 
GSL_VAR const char * ool_version;

__END_DECLS

#endif /*__OOL_VERSION_H__*/

/*----------------------------------------------------------------------------*/

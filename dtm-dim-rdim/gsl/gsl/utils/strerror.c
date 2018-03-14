/*
 * strerror.c --- ANSI C compatible system error routine
 */

/* 
 * Copyright (C) 1986, 1988, 1989, 1991 the Free Software Foundation, Inc.
 * From gawk.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 */

#if 0
#include <stdio.h>
#endif

extern int sys_nerr;
extern char *sys_errlist[];

char *
strerror(n)
int n;
{
        static char mesg[30];

        if (n < 0 || n >= sys_nerr) {
                sprintf(mesg, "Unknown error (%d)", n);
                return mesg;
        } else
                return sys_errlist[n];
}

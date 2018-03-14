/* xmalloc.c -- safe versions of malloc and realloc */

/* This file is part of GNU Info, a program for reading online documentation
   stored in Info format.

   This file has appeared in prior works by the Free Software Foundation;
   thus it carries copyright dates from 1988 through 1993.

   Copyright (C) 1988, 1989, 1990, 1991, 1992, 1993 Free Software
   Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

   Written by Brian Fox (bfox@ai.mit.edu). */

#if !defined (ALREADY_HAVE_XMALLOC)
#include <stdio.h>
#include <sys/types.h>

extern void *malloc (int), *realloc (void *, int);
static void memory_error_and_abort (const char *);

extern void *xmalloc (int), *xrealloc (void *, int);


/* **************************************************************** */
/*                                                                  */
/*                 Memory Allocation and Deallocation.              */
/*                                                                  */
/* **************************************************************** */

/* Return a pointer to free()able block of memory large enough
   to hold BYTES number of bytes.  If the memory cannot be allocated,
   print an error message and abort. */

void *
xmalloc (int bytes)
{
  void *temp = malloc (bytes);

  if (!temp)
    memory_error_and_abort ("xmalloc");
  return (temp);
}

void *
xrealloc (void *pointer, int bytes)
{
  void *temp;

  if (!pointer)
    temp = malloc (bytes);
  else
    temp = realloc (pointer, bytes);

  if (!temp)
    memory_error_and_abort ("xrealloc");

  return (temp);
}

static void
memory_error_and_abort (const char * fname)
{
  fprintf (stderr, "%s: Out of virtual memory!\n", fname);
  abort ();
}
#endif /* !ALREADY_HAVE_XMALLOC */

#! /bin/sh

# SPDX-FileCopyrightText: 2019 - 2023 sudorook <daemon@nullcodon.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

set -e

test -z "$srcdir" && srcdir=$(dirname "$0")
test -z "$srcdir" && srcdir=.

cd "$srcdir"

autoreconf --force --install --warnings=all

if test -z "${NOCONFIGURE}"; then
  set -x
  "$srcdir"/configure "$@"
fi

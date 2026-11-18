#!/bin/bash
# setRetroMoCoBoxHome.sh
#
# Sourced by the cluster_submit*.sh scripts to set RETROMOCOBOX_HOME to the
# installed copy that matches the version defined in retroMocoBoxVersion.m
# (in the repo root, one level above this 'cluster' folder).
#
# e.g. out = '1.0.3dev';  ->  /cubric/software/matlab.toolboxes/retroMoCoBox_v1.0.3
#
# code generated for D. Gallichan by Claude, date: 18/11/26

RETROMOCOBOX_INSTALL_ROOT=/cubric/software/matlab.toolboxes

# folder containing this file (follows symlinks, works when sourced)
_rmcb_clusterDir="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
_rmcb_versionFile="${_rmcb_clusterDir}/../retroMocoBoxVersion.m"

if [ ! -f "${_rmcb_versionFile}" ]; then
    echo "ERROR: could not find ${_rmcb_versionFile}" >&2
    return 1
fi

# pull the string out of the line:  out = '1.0.3dev';
RETROMOCOBOX_VERSION_FULL=$(sed -n "s/^[[:space:]]*out[[:space:]]*=[[:space:]]*'\([^']*\)'.*/\1/p" "${_rmcb_versionFile}" | head -n 1)

if [ -z "${RETROMOCOBOX_VERSION_FULL}" ]; then
    echo "ERROR: could not read version string from ${_rmcb_versionFile}" >&2
    return 1
fi

# installed folders are named by the numeric part only (drop e.g. 'dev')
RETROMOCOBOX_VERSION=$(echo "${RETROMOCOBOX_VERSION_FULL}" | grep -oE '^[0-9]+(\.[0-9]+)*')

export RETROMOCOBOX_HOME=${RETROMOCOBOX_INSTALL_ROOT}/retroMoCoBox_v${RETROMOCOBOX_VERSION}

if [ ! -d "${RETROMOCOBOX_HOME}" ]; then
    echo "ERROR: retroMoCoBox version ${RETROMOCOBOX_VERSION_FULL} (from ${_rmcb_versionFile})" >&2
    echo "       expected to be installed at ${RETROMOCOBOX_HOME}, but that folder does not exist." >&2
    echo "       Available versions:" >&2
    ls -d ${RETROMOCOBOX_INSTALL_ROOT}/retroMoCoBox_v* 2>/dev/null | sed 's/^/         /' >&2
    return 1
fi

echo "Using retroMoCoBox version ${RETROMOCOBOX_VERSION_FULL}: ${RETROMOCOBOX_HOME}"

unset _rmcb_clusterDir _rmcb_versionFile

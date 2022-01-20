#!/usr/bin/env bash

# Jupyter is unable to apply the syntax highlighting colors from
# custom themes in jupyterthemes.
# By saving the notebook in the browser with save as the colors
# are preserved but the header remains and the scrolling will be
# disabled.
# These issues can be resolved the way below.

if [ ! -f "$1" ];
then
	echo "Usage: $0 <notebook_saved_from_browser.html>"
	exit -1
fi

sed -iE ':a;s/<!--.*-->//g;/<!--/{N;ba};'\
's/@font-face {.*}//g;'\
's/<html/<html style="height:100%;"/;'\
's/<body/<body style="height:100%;overflow-y:auto!important;"/;'\
's/\(<div id="\?header"\?\)/<!-- \1/;'\
's/\(<div id="\?site"\? style="\?.*"\?>\)/\1 -->/;'\
's/\(id="?notebook_panel"?\) style="overflow:auto!important;"/\1/g;'\
's/\(id="?notebook_panel"?\)/\1 style="overflow:auto!important;"/g;'\
 "$1"

#!/bin/bash
find . -iname "log" | xargs -I {} grep '^INSERT' {} > to_db.sqlite
find $EXPERIMENTS -iname "*log" | xargs grep -L '^INSERT' | xargs -I {} tail -n1 -v {} > to_be_run_again.txt


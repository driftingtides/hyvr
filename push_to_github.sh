#!/bin/bash

# make sure we are in the mercurial repo
if ! hg status 1> /dev/null 2>&1
then
    echo "You are not in a mercurial repository. Make sure to start the script from the top level of
    the mercurial repo."
    exit 1
fi

# make sure we are on the top level 
if [ ! -f .hgignore ]
then
    echo "Make sure you are at the top level of the repo."
    exit 1
fi

# optionally add a tag
while getopts "h?t:" opt; 
do
    case "$opt" in 
        h|\?)
            echo "Usage: push_to_github.sh [-t tagname]"
            echo ""
            echo "This script pushes the current content of the mercurial repo to github."
            echo "This should only be run from the top level of the mercurial repo."
            echo "If you want to add a tag, use the -t option and supply a tagname"
            exit 0
            ;;
        t)
            tag=1
            tagname=$OPTARG
            ;;
    esac
done


# checkout repo in temporary directory
cd ..
mkdir tmp_git_repo
cd tmp_git_repo
git clone git@github.com:driftingtides/hyvr.git

# add/remove files
cd hyvr
rm -rf * .gitignore
cp -r ../../hyvr/* ../../hyvr/.hgignore .
mv .hgignore .gitignore
echo "push_to_github.sh" >> .gitignore
echo ".hgtags" >> .gitignore
echo "pytest.ini" >> .gitignore
git add -A

# commit
git commit
if [ "$tag" ]
then
    git tag "$tagname"
fi
# push to github
git push origin "$tagname"

# push to github pages
cd docs
make gh-pages


# cleanup
cd ../../../
rm -rf tmp_git_repo

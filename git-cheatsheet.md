# Basics

To create a new repository on your machine, make sure your computer has an ssh key that has been [uploaded to github](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account), and then

    git clone git@github.com:roed314/padic_square.git

This will create a `padic_square` folder in your current directory, and all future commands should be executed within that folder.

To pull in changes from github that others have made, do

    git pull

To contribute changes, there are two steps.  First, you have to commit your changes locally so that git knows about them (Message should be a short description of what you're doing).

    git commit -a -m "Message"

This only works when you're editing files that git already knows about.  If you want to add new files, instead you should do

    git add -A
    git commit -m "Message"

After this, you can push your changes to github

    git push

Note that this will fail if someone else has pushed since you last pulled.  To proceed, just run `git pull` and then `git push`.  If you're editing the same file as someone else, there may be conflicts.  Open the files where git tells you that there are conflicts, edit them to resolve the conflicts, do `git commit -m "Merging"` and then `git push`.  Feel free to seek help in this case.

# Other commands

You can see the current status of git using

    git status

If you're editing a file, you can see the changes that have been made since the last commit using

    git diff

If there are files that you want git to ignore, you can edit the `.gitignore` file in this repository.

Git has lots more features.  Here's a slightly longer [cheatsheet](https://www.atlassian.com/git/tutorials/atlassian-git-cheatsheet) and a [tutorial](https://www.atlassian.com/git).
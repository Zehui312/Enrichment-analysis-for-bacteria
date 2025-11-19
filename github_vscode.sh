#=================================================================
#+++++++++++++++++++++++Step How to get github author ++++++++++++
#=================================================================
ssh-keygen -t ed25519 -C "zyu@stjude.org"

eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519

cat ~/.ssh/id_ed25519.pub # add this key at GitHub → Settings → SSH and GPG keys


git remote set-url origin git@github.com:Zehui312/Enrichment-analysis-for-bacteria.git
ssh -T git@github.com

# git clone git@github.com:Zehui312/BacDrop-analysis-Pipeline.git
git clone https://github.com/Zehui312/Enrichment-analysis-for-bacteria.git
#=================================================================
#+++++++++++++++++++++++Step update github +++++++++++++++++++++++
#=================================================================
ssh -T git@github.com
git push origin main
git remote set-url origin git@github.com:Zehui312/Enrichment-analysis-for-bacteria.git
git push origin main

git commit -am "update tracked files"
git push origin main

git add .
git commit -m "update project files"
git push origin main

#test
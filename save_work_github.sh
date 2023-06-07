#!/bin/bash


echo Here is your status

git status

git add .

echo Work added

read comment

echo $comment

git commit -m "$comment"

echo Work committed

git push origin main

echo Work pushed

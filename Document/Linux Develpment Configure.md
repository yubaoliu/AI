---
Title: Linux Develpment Configure
Date: Nov 1 2018
Author: yubaoliu89@gmail.com
---

# Develpment Configure

## CPP

```cpp
sudo apt-get install build-essential
```
#Package Management

## List Packges

```sh
dpkg --list
dpkg --list | grep <package_name>
```

## Uninstall packages

```sh
sudo apt-get --purge remove <program>
```



## Uninstall Default Program

1. uninstall LibreOffice

    ```sh
    sudo apt-get remove --purge libreoffice*
    sudo apt-get clean
    sudo apt-get autoremove
    ```

1.  Uninstall thunderbird

    ```sh
    sudo apt-get remove --purge thunderbird*
    ```



# Ubuntu Issues

## Software Center

[Fix Ubuntu Software Center not loading issue in Ubuntu 16.04 LTS](https://www.fosslinux.com/1768/fix-ubuntu-software-center-not-loading-issue-in-ubuntu-16-04-lts.htm/)

Fix Ubuntu 16.04 Software Center not loading apps issue:

Step 1) Launch ‘Terminal’.

Step 2) Run the following command to update the repository sources.

```
sudo apt-get update
```

Step 3) Now install the updates.

```
sudo apt-get upgrade
```


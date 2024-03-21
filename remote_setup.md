# Things to do for LCC3

## Generate an SSH-Key for passwordless login

1. Generate a key pair: `ssh-keygen -t ed25519 -C "some description (e.g. e-mail address, hostname, etc.)"`
2. Copy the public key to the destination host: `ssh-copy-id -i ~/.ssh/id_ed25519.pub <username>@lcc3.uibk.ac.at`

- Source: <https://www.ssh.com/academy/ssh/copy-id#setting-up-public-key-authentication>

## Setup Visual Studio Code for Remote SSH work

### Configure SSH Shortcut

1. Open the SSH configuration
    1. Add a new entry like this:

    ````ssh-config
    Host lcc3
      HostName lcc3.uibk.ac.at
      User <username>
      IdentityFile ~/.ssh/id_ed25519
    ````

### Configure LCC3-specific settings in Remote SSH extension

1. Open Visual Studio Code
2. Install the Remote SSH extension
3. Open the Remote SSH settings
    1. Add a new entry in `Remote.SSH: Server Install Path`
        - Key: `lcc3`
        - Value: `/scratch/<username>`

### Configure clangd stuff

1. Install the clangd VS Code extension on LCC3
2. When asked, have VS Code download and install clangd
3. 

cmake stuff

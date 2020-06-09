# Seng440 project

A project on Singular Value Decomposition (SVD) and optimizing its performance on an arm machine.

## Requirements

- qemu
- qemu-kvm
- qemu-system-arm
- libvirt-clients
- libvirt-daemon-system
- bridge-utils
- virt-manager
  - arch: arm
  - machine type: virt-2.11
  - Fedora 29
  - 2GB Ram
  - 1 CPU
- gcc compiler

## Installation

Create a VM in virt-manager using the settings specified in the [requirements](#requirements) section.

Add the following kernel arguments:

```bash
console=ttyAMA0
rw
root=LABEL=_/
rootwait
ipv6.disable=1
```

## Documentation

To see the documentation go to the `docs` folder.

The implementation is based on a set of [slides provided by the professor](docs/Embedded_Systems_Slides_WRAPON_lesson_112.pdf) in the `docs` folder.

## Authors

- Robert Tulip - rtulip@uvic.ca
- Michail Roesli - mroesli@uvic.ca

CREATE TABLE `perAtomASA` (
`perAtomASA_id` INTEGER NOT NULL ,
`PDB` TEXT NOT NULL DEFAULT '',
`Chain` TEXT NOT NULL DEFAULT '',
`Residue` TEXT NOT NULL DEFAULT '',
`ResId` INTEGER NOT NULL,
`iCode` TEXT NOT NULL DEFAULT '',
`Symbol` TEXT NOT NULL DEFAULT '',
`Atom` TEXT NOT NULL DEFAULT '',
`ASA` float NOT NULL,
`Bfactor` float NOT NULL,
`Seperated` char(1) NOT NULL DEFAULT '0',
PRIMARY KEY (`perAtomASA_id`)
);
DROP TABLE IF EXISTS `Ndrieding`;

CREATE TABLE `Ndrieding` (
  `fullDrieding_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(4) NOT NULL DEFAULT '',
  `DonorChain` varchar(2) NOT NULL DEFAULT '',
  `DonorResId` int(11) NOT NULL,
  `DonorICode` varchar(2) NOT NULL DEFAULT '',
  `DonorSymbol` varchar(4) NOT NULL DEFAULT '',
  `AccChain` varchar(2) NOT NULL DEFAULT '',
  `AccResId` int(11) NOT NULL,
  `AccICode` varchar(2) NOT NULL DEFAULT '',
  `AccSymbol` varchar(4) NOT NULL DEFAULT '',
  `Energy` float NOT NULL,
  PRIMARY KEY (`fullDrieding_id`)
)

DROP TABLE IF EXISTS `interfaceDist`;
CREATE TABLE `interfaceDist` (
  `interfaceDist_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(5) NOT NULL DEFAULT '',
  `Chains` varchar(6) NOT NULL DEFAULT '',
  `Chain` varchar(4) NOT NULL DEFAULT '',
  `ResId` int(11) NOT NULL,
  `iCode` varchar(2) NOT NULL DEFAULT '',
  `Symbol` varchar(4) NOT NULL DEFAULT '',
  `Atom` varchar(2) NOT NULL DEFAULT '',
  `MinDist` float NOT NULL,
  PRIMARY KEY (`interfaceDist_id`)
)

DROP TABLE IF EXISTS `NinterfaceAtoms`;

CREATE TABLE `NinterfaceAtoms` (
  `PDB` varchar(5) NOT NULL DEFAULT '',
  `Chain` varchar(5) NOT NULL DEFAULT '',
  `Residue` varchar(4) NOT NULL DEFAULT '',
  `ResId` int(11) NOT NULL,
  `iCode` varchar(2) NOT NULL DEFAULT '',
  `Symbol` varchar(4) NOT NULL DEFAULT '',
  `atom` varchar(2) NOT NULL DEFAULT '',
  `diffASA` double NOT NULL DEFAULT '0',
  `pk` int(11) NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`pk`)
)

DROP TABLE IF EXISTS `electrostat`;
CREATE TABLE `electrostat` (
  `electrostat2_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(4) NOT NULL DEFAULT '',
  `electro` float NOT NULL,
  `pp` int NOT NULL,
  `mm` int NOT NULL,
  `pm` int NOT NULL,
  PRIMARY KEY (`electrostat2_id`)
)
DROP TABLE IF EXISTS `nRMSD`;
CREATE TABLE `nRMSD` (
  `nRMSD_id` int(11) NOT NULL AUTO_INCREMENT,
  `Complex` varchar(5) NOT NULL DEFAULT '',
  `Unbound` varchar(5) NOT NULL DEFAULT '',
  `Chain` varchar(5) NOT NULL DEFAULT '',
  `UnboundChain` varchar(5) NOT NULL DEFAULT '',
  `RMSD` float NOT NULL,
  `iRMSD` float NOT NULL,
  `Atoms` int(11) NOT NULL,
  `iAtoms` int(11) NOT NULL,
  PRIMARY KEY (`nRMSD_id`)
)
DROP TABLE IF EXISTS `interfaceVDW`;
CREATE TABLE `interfaceVDW` (
  `interfaceVDW_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(5) NOT NULL DEFAULT '',
  `VDV` float NOT NULL,
  `VDVx6` float NOT NULL,
  `ClashV` float NOT NULL,
  `ClashS` float NOT NULL,
  PRIMARY KEY (`interfaceVDW_id`)
) 
DROP TABLE IF EXISTS `interfacePeriphrial`;
CREATE TABLE `interfacePeriphrial` (
  `interfacePeriphrial_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(5) NOT NULL DEFAULT '',
  `Chain` varchar(5) NOT NULL DEFAULT '',
  `ResId` int(11) NOT NULL,
  `Symbol` varchar(5) NOT NULL DEFAULT '',
  `Peri` float NOT NULL,
  `PropPeri` float NOT NULL,
  PRIMARY KEY (`interfacePeriphrial_id`)
) 
DROP TABLE IF EXISTS `proteinComplex`;
CREATE TABLE `proteinComplex` (
  `proteinComplex_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(5) NOT NULL DEFAULT '',
  `UnboundChainA` varchar(6) NOT NULL DEFAULT '',
  `NameA` varchar(55) NOT NULL DEFAULT '',
  `UnboundChainB` varchar(6) NOT NULL DEFAULT '',
  `NameB` varchar(55) NOT NULL DEFAULT '',
  PRIMARY KEY (`proteinComplex_id`)
)


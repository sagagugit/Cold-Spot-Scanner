DROP TABLE IF EXISTS `donors2`;
CREATE TABLE `donors2` (
  `donors2_id` int(11) NOT NULL AUTO_INCREMENT,
  `Residue` varchar(4) NOT NULL DEFAULT '',
  `Symbol` varchar(5) NOT NULL DEFAULT '',
  `Hydrogen` varchar(5) NOT NULL,
  PRIMARY KEY (`donors2_id`)
) 
LOCK TABLES `donors2` WRITE;
INSERT INTO `donors2` VALUES (18,'-','N','H'),(13,'ARG','NE','HE'),(9,'ARG','NH1','1HH1'),(10,'ARG','NH1','2HH1'),(11,'ARG','NH2','1HH2'),(12,'ARG','NH2','2HH2'),(14,'ASN','ND2','1HD2'),(15,'ASN','ND2','2HD2'),(7,'CYS','SG','HG'),(17,'GLN','NE2','1HE2	'),(8,'HIS','NE2','HE2'),(4,'LYS','NZ','1HZ'),(5,'LYS','NZ','2HZ'),(6,'LYS','NZ','3HZ'),(2,'SER','OG','HG'),(3,'THR','OG1','HG1'),(16,'TRP','NE1','HE1'),(1,'TYR','OH','HH');
UNLOCK TABLES;

function toFile(sc, fname)
%converts segmentation cluster or group of segmentation clusters to 
%a YAML string representation which it saves to disk
%function toFile(sc, fname)


import('org.yaml.snakeyaml.Yaml');
yaml = Yaml();
outfile = java.io.FileWriter(fname);
yaml.dump(sc.toHashMap, outfile);
outfile.close();
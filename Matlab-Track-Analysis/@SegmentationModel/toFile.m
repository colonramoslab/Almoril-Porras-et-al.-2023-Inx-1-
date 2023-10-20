function toFile(sm, fname)
%converts segmentation model or group of segmentation models to 
%a YAML string representation which it saves to disk
%function toFile(sc, fname)


import('org.yaml.snakeyaml.Yaml');
yaml = Yaml();
outfile = java.io.FileWriter(fname);
yaml.dump(sm.toHashMap, outfile);
outfile.close();
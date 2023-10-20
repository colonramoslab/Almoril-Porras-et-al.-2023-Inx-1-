function sc = fromFile(fname, varargin)
%converts from a string representation to the array of segmentation
%clusters
%function sc = fromString(str, eset)
import('org.yaml.snakeyaml.Yaml');
yaml = Yaml();
is = java.io.FileInputStream(fname);
obj = yaml.load(is);
sc = SegmentationCluster.fromHashMap(obj, varargin{:});
is.close();
function sm = fromString(s, varargin)
%converts from a string representation to a segmentation model
%function sm = fromString(s, varargin)
import('org.yaml.snakeyaml.Yaml');
yaml = Yaml();
obj = yaml.load(s);
sm = SegmentationModel.fromHashMap(obj, varargin{:});

function filename = test_phasesym(n)
  %n = number of iterations
  
  %------------------------------------------------------%
  % Input parameter initialization
  %------------------------------------------------------%

  % Image array initialization
  im_init = cell(1,5);
  im_init{1} = 'data/cameraman.tif';
  im_init{2} = 'data/rice.png';
  im_init{3} = 'data/coins.png';
  im_init{4} = 'data/glass.png';
  im_init{5} = 'data/blobs.png';

  % No of Scales
  scale_init = [3:1:6];

  % No of filter Orientations
  orient_init = [4:1:6];

  % Wavelength of smallest scale filter
  minWaveLength_init = 3

  % scaling factor between successive filters
  mult_init = 2.1

  % ratio of std deviation of gaussian for log Gabors
  sigmaOnf_init = 0.55

  % Ratio of angular interva between filter orientations
  dThetaOnSigma_int = 1.2

  % no of std deviations of noise energy beyond the mean
  k_init = [1:1:3];

  % polarity
  polarity_init = [-1:1:1];

  for i = 1:n % for total number of interations
    im = im_init{ceil(1 + (5-1).*rand(1,1))};
    image = imread(im);
    image = grayscale(image);
    scale = scale_init(ceil(0+(4-0).*rand(1,1)));
    orient = orient_init( );
    minWaveLength = minWaveLength_init;
    mult = mult_init;
    sigmaOnf = sigmaOnf_init;
    dThetaOnSigma = dThetaOnSigma_int;
    k = k_init(ceil(0+(3-0).*rand(1,1)));
    polarity = polarity_init(ceil(0+(3-0).*rand(1,1)));

    # use proper filename here
    filename = strcat('file',int2str(i), im, int2str(scale), int2str(orient),
			    int2str(minWaveLength), int2str(mult), int2str(sigmaOnf),
			    int2str(dThetaOnSigma), int2str(k), int2str(polarity));

    
    [phasesym, orientation, totalEnergy] = phasesym(im, scale, orient, minWaveLength, mult,
						    sigmaOnf, dThetaOnSigma, k, polarity);

    save(filename_input, 'im', 'scale', 'orient', 'minWaveLength', 
	 'mult', 'sigmaOnf', 'dThetaOnSigma', 'k', 'polarity',
	 'phasesym', 'orientation', 'totalEnergy');

  end
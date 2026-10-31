#ifndef GR_SAMPLINGENTRYBUILDER_H
#define GR_SAMPLINGENTRYBUILDER_H

class SamplingEntryBuilder {

  SamplingEntryBuilder();

  SamplingEntryBuilder(const SamplingEntryBuilder&);
  SamplingEntryBuilder&
  operator=(const SamplingEntryBuilder&);

public:

  virtual
  ~SamplingEntryBuilder();

  static SamplingEntryBuilder*
  instance();

};

#endif

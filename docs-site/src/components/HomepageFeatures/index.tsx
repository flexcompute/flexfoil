import type {ReactNode} from 'react';
import clsx from 'clsx';
import Heading from '@theme/Heading';
import styles from './styles.module.css';

type FeatureItem = {
  title: string;
  description: ReactNode;
};

const FeatureList: FeatureItem[] = [
  {
    title: 'Python API',
    description: (
      <>
        <code>pip install flexfoil</code> — solve airfoils in three lines of
        code. Run polar sweeps, export to pandas, and plot with matplotlib.
      </>
    ),
  },
  {
    title: 'Web App',
    description: (
      <>
        Interactive browser-based analysis with flow visualization, polar
        plots, a data explorer, plot builder, and inverse design tools.
      </>
    ),
  },
  {
    title: 'XFOIL Heritage',
    description: (
      <>
        Faithful reimplementation of Mark Drela&apos;s XFOIL solver —
        viscous-inviscid coupling, e<sup>N</sup> transition, and panel method.
      </>
    ),
  },
];

function Feature({title, description}: FeatureItem) {
  return (
    <div className={clsx('col col--4')}>
      <div className="text--center padding-horiz--md">
        <Heading as="h3">{title}</Heading>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures(): ReactNode {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}

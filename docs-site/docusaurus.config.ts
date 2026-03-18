import {themes as prismThemes} from 'prism-react-renderer';
import type {Config} from '@docusaurus/types';
import type * as Preset from '@docusaurus/preset-classic';

const config: Config = {
  title: 'FlexFoil',
  tagline: 'Airfoil analysis powered by XFOIL — in Python and the browser',
  favicon: 'img/logo.svg',

  future: {
    v4: true,
  },

  url: 'https://foil.flexcompute.com',
  baseUrl: '/docs/',

  organizationName: 'flexfoil',
  projectName: 'flexfoil',

  onBrokenLinks: 'throw',

  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  // Enable Mermaid for diagrams
  markdown: {
    mermaid: true,
  },
  themes: ['@docusaurus/theme-mermaid'],

  presets: [
    [
      'classic',
      {
        docs: {
          routeBasePath: '/',
          sidebarPath: './sidebars.ts',
          editUrl: 'https://github.com/flexfoil/flexfoil/tree/main/docs-site/',
        },
        blog: false,
        theme: {
          customCss: './src/css/custom.css',
        },
      } satisfies Preset.Options,
    ],
  ],

  themeConfig: {
    // image: 'img/social-card.jpg',
    colorMode: {
      respectPrefersColorScheme: true,
    },
    mermaid: {
      theme: {light: 'neutral', dark: 'dark'},
    },
    navbar: {
      title: 'FlexFoil',
      logo: {
        alt: 'FlexFoil Logo',
        src: 'img/logo.svg',
      },
      items: [
        {
          type: 'docSidebar',
          sidebarId: 'tutorialSidebar',
          position: 'left',
          label: 'Documentation',
        },
        {
          href: 'https://github.com/flexfoil/flexfoil',
          label: 'GitHub',
          position: 'right',
        },
      ],
    },
    footer: {
      style: 'dark',
      links: [
        {
          title: 'Docs',
          items: [
            {
              label: 'Get Started',
              to: '/intro',
            },
            {
              label: 'Python API',
              to: '/python-api',
            },
            {
              label: 'Web App',
              to: '/web-app',
            },
          ],
        },
        {
          title: 'More',
          items: [
            {
              label: 'Solver Internals',
              to: '/internals/xfoil-architecture-spec',
            },
            {
              label: 'GitHub',
              href: 'https://github.com/flexfoil/flexfoil',
            },
          ],
        },
      ],
      copyright: `Copyright © ${new Date().getFullYear()} FlexFoil Project. Built with Docusaurus.`,
    },
    prism: {
      theme: prismThemes.github,
      darkTheme: prismThemes.dracula,
      additionalLanguages: ['fortran', 'rust', 'toml'],
    },
  } satisfies Preset.ThemeConfig,
};

export default config;
